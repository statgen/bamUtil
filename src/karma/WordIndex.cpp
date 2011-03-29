/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "InplaceMerge.h"
#include "Performance.h"
#include "WordIndex.h"
#include "Error.h"
#include "Util.h"

#include <ctype.h>

#include <iomanip>
#include <iostream>
#include <map>

//
// The WordIndex class exists to rapidly map an integer key to a list
// of genome positions.
//
// The integer key is simply the bitwise concatenation of the base pair
// index values (0-3), taken two bits at a time.
//
// So, the bit value looks like this:
//
// 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00
//       ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//  0  0  0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3   0-3
// For example:
// ATCGATCGATCGATC will look like this:
// base  4: 012301230123012
// base  2: 00 0110 1100 0110 1100 0110 1100 0110
// base 10: 113690310
// base 16: 06C6C6C6
//
//

String WordIndex::IntegerToBPSeq(wordInteger_t word)
{
    String sequence("");

    if (word > wordsCountIndexMask)
        error("%d needs to be a non-negative integer < wordsCountIndexMask\n", word);

    for (unsigned int i = 0; i < wordSize; i ++)
        sequence += "N";

    for (unsigned int i = 0; i < wordSize; i++)
    {
        sequence[wordSize-1-i] = GenomeSequence::int2base[word & 3];
        word >>= 2;
    }
    return sequence;
}


WordIndex::WordIndex()
{
    numberWords = 0;
    numberWordPositions = 0;
    hashindices = NULL;
    wordCounts = NULL;
    wordReachedCutoffBitvector = NULL;
    wordReachedCutoffBitvectorLength = 0;

    debugFlag = false;
    verboseFlag = false;
    occurrenceCutoff = 0;
}

WordIndex::~WordIndex()
{
    hashIndicesMemoryMap.close();
    wordPositionsMemoryMap.close();
}

bool WordIndex::AllocateWordsCountHashMemory()
{
    assert(wordSize <= 16 && wordSize > 0);

    // 5 -> 4 for the header, 1 for an extra word at the end of hashindices.
    if (hashIndicesMemoryMap.create(hashIndicesFilename, sizeof(unsigned int) *(wordsCountHashSize + 5)))
    {
        return false;
    }

    *((unsigned int *) hashIndicesMemoryMap.data + 0) = UMWI_HI_COOKIE;
    *((unsigned int *) hashIndicesMemoryMap.data + 1) = UMWI_VERSION;
    *((unsigned int *) hashIndicesMemoryMap.data + 2) = wordSize;
    //  *((uint64_t *) hashIndicesMemoryMap.data + 3) = wordsCountHashSize;
    hashindices = (unsigned int *) hashIndicesMemoryMap.data + 3;

    // wordCounts is only needed during ::create().  It is
    // deleted in ::StoreWordPositions().
    wordCounts = new unsigned int [wordsCountHashSize];
    if (wordCounts == NULL) return false;

    for (uint64_t i = 0; i < wordsCountHashSize; i ++)
    {
        wordCounts[i] = 0;
    }
    return true;
}

std::ostream &operator << (std::ostream &stream, WordIndex &wordIndex)
{
    stream << "File: " << wordIndex.hashIndicesFilename << "\n";
    stream << "wordSize: " << wordIndex.wordSize << "\n";
    stream << "wordsCountHashSize: " << wordIndex.wordsCountHashSize << "\n";
    stream << "wordsCountIndexMask: " << wordIndex.wordsCountIndexMask << "\n";
    stream << "numberWords: " << wordIndex.numberWords << "\n";
    stream << "numberWordPositions: " << wordIndex.numberWordPositions << "\n";
    return stream;
}



void WordIndex::WordsCountHash(GenomeSequence & gs)
{
    assert(wordSize <= 16 && wordSize > 0);

    std::cout<< "Finding all words of size " << wordSize << std::endl << std::flush;

    if (!AllocateWordsCountHashMemory())
        error("Memory allocation for words count failed\n");

    numberWords = 0; // including repeats;

    wordInteger_t sequenceInteger = 0;
    unsigned int bpCounter = 0;
    unsigned int nextprogress = 10000000;

    for (unsigned int basePairIndex = 0; basePairIndex < gs.sequenceLength(); basePairIndex ++)
    {
        // XXX refactor this line:
        uint8_t addBaseInteger = gs.getInteger(basePairIndex);
        //
        // For the time being, we drop all words that have unknown encodings.
        // These are 'N' or 'M' values in the reference Genome.
        // 'N' ==> iNdelibles
        // 'M' ==> the reference has some of these, but I don't know what they mean
        //
        if (addBaseInteger > 3)
        {
            // recounting up to <wordSize> bp;
            bpCounter = 0;
            sequenceInteger = 0;
            continue;
        }

        sequenceInteger = ((sequenceInteger<<2) |  addBaseInteger) & wordsCountIndexMask;

        if (++bpCounter >= wordSize)
        {
            if (basePairIndex >= nextprogress)
            {
                int megabase = basePairIndex / 1000000;
                std::cout << megabase << " Mb\r" << std::flush;
                nextprogress += 10000000;
            }

            if (wordCounts[sequenceInteger]==0) numberWords ++ ;

            if (wordCounts[sequenceInteger] < WORDCOUNT_MAX)
            {
                wordCounts[sequenceInteger] ++;
            }
        } //

    } // end of for unsigned int basePairIndex

    std::cout << "\nTotal " << numberWords << " words of size " << wordSize << std::endl << std::flush;

}

//
// XXX for debug checking
//
void WordIndex::WriteWordsCountHash(const String & prefix)
{
#if 0
    if (wordSize > 11)
    {
        FILE * highrepeatwordcounts = fopen(prefix + ".highrepeat.counts", "wt");
        FILE * lowrepeatwordcounts = fopen(prefix + ".lowrepeat.counts", "wt");
        if (highrepeatwordcounts == NULL || lowrepeatwordcounts == NULL)
            error("Failed to open output file for word counts\n");

        fprintf(highrepeatwordcounts, "word count\n");
        fprintf(lowrepeatwordcounts, "word count\n");
        for (unsigned int i = 0; i < wordsCountHashSize; i ++)
        {
            if (hashindices[i] != UINT32_MAX)
            {
                String word = IntegerToSeq(hashindices[i]);

                if (wordCounts[i] > occurrenceCutoff)
                {
                    word.Write(highrepeatwordcounts);
                    fprintf(highrepeatwordcounts, " %d\n", wordCounts[i]);
                }
                else
                {
                    word.Write(lowrepeatwordcounts);
                    fprintf(lowrepeatwordcounts, " %d\n", wordCounts[i]);
                }
            }
        }
    }
#endif
}


//
// Initialize hashindices.
// Each array entry will now have a value assigned as follows:
// <integer>    -> an integer index into the WordPositions array.
//
// hashindices[i+1] - hashindices[i] ==> number of occurrences of this word.
//
void WordIndex::IndexPositionsForLowRepeatWords()
{
    uint64_t    wordsExceedingMaxRepeatCount = 0;
    uint64_t    genomeLocationsOmittedCount = 0;
    assert(wordSize < 16 && wordSize > 0);

    unsigned int numberWordsIndexed = 0;

    // wordReachedCutoffBitvectorLength gets set in ::setWordSize()
    wordReachedCutoffBitvector =
        (unsigned char *) malloc(wordReachedCutoffBitvectorLength);

    memset(wordReachedCutoffBitvector, 0, wordReachedCutoffBitvectorLength);

    numberWordPositions = 0;
    for (uint64_t i = 0; i < wordsCountHashSize; i ++)
    {
        hashindices[i] = numberWordPositions;
        // the word does not show up in the reference genome:
        if (wordCounts[i]==0) continue;

        // word occurs too often:
        if (wordCounts[i] > occurrenceCutoff)
        {
            genomeLocationsOmittedCount += wordCounts[i];
            wordsExceedingMaxRepeatCount++;
            // mark it as having reached cutoff
            wordReachedCutoffBitvector[i>>3] |= (1<<(i&0x7));
            continue;
        }

        // word positions have to include their count, so add one here.
        // We need the count so we can later sort the set of word positions.
        numberWordPositions += wordCounts[i];
        numberWordsIndexed ++;
    }
    hashindices[wordsCountHashSize] = hashindices[wordsCountHashSize - 1];

    std::cout << wordsExceedingMaxRepeatCount << " words were omitted due to exceeding " << occurrenceCutoff << " cutoff." << std::endl;
    std::cout << genomeLocationsOmittedCount << " genome match positions were omitted due to those words exceeding the cutoff." << std::endl;
    std::cout << "These positions will be stored in a secondary hash." << std::endl;

    std::cout << "\nTotal " << numberWordPositions << " positions for " << numberWordsIndexed
              << " repeats occuring <= " <<  occurrenceCutoff << " times" << std::endl;
}

// to be completed
int WordIndex::GetOccurrenceCutoff(double pctWordsToKeep)
{
    assert(wordSize <= 16 && wordSize > 0);
    assert(pctWordsToKeep>0.0);

    //  unsigned int indexFromHighRepeats = (unsigned int) numberWords * (1.0 - pctWordsToKeep);
    int cutoff = 1;

    // sort XXWordCounts to get the cutoff BUT do NOT distort the original order;
    //QuickSort();
    return cutoff;

}

void WordIndex::IndexPositionsForLowRepeatWords(double pctWordsToKeep)
{
    int cutoff = GetOccurrenceCutoff(pctWordsToKeep);
    occurrenceCutoff = cutoff;
    IndexPositionsForLowRepeatWords();

}

void WordIndex::StoreWordPositions(GenomeSequence & gs, WordHash &wordHashLeft, WordHash &wordHashRight)
{

    assert(wordSize <= 16);

    assert(numberWordPositions > 0);

    //
    // we could memory map a file, then write into it, but because the writes
    // to the wordpositions array are essentially random, and because it is so
    // large, the filesystem gets bogged down, the blocks while trying to write
    // dirty pages.  So mmap() does not help in this instance.  Use malloc()
    // to simply get normal dynamically allocated memory... we'll do normal file I/O at
    // the end of this routine.
    //
    wordpositions = (unsigned int *) malloc(sizeof(unsigned int) * numberWordPositions);

    // all values should get set, but this is safer
    for (unsigned int i = 0; i < numberWordPositions; i ++)
        wordpositions[i] = UINT32_MAX;

    uint64_t sequenceInteger = 0;

    unsigned int bpCounter = 0;
    unsigned int nextprogress = 10000000;

    std::cout << "Storing positions for non-highly repeated words - table size = "
              << (float)numberWordPositions * sizeof(int) / 1024.0 / 1024.0 / 1024.0
              << " Gbytes" << std::endl;

    // again, we loop over all words in the genome.  For each
    // word found here, we need to add it to the list of words
    // in wordpositions.
    for (unsigned int basePairIndex = 0; basePairIndex < gs.sequenceLength(); basePairIndex ++)
    {
        uint8_t addBaseInteger = gs.getInteger(basePairIndex);

        // if an invalid basepair, start over from scratch
        if (addBaseInteger > 3)
        {
            // recounting up to <wordSize> bp;
            bpCounter = 0;
            sequenceInteger = 0;
            continue;
        }

        sequenceInteger = ((sequenceInteger<<2) |  addBaseInteger) & wordsCountIndexMask;

        // when we have collected "wordSize" basepairs into our integer, we'll
        // add the position for this word into the wordpositions table.
        //
        if (++bpCounter >= wordSize)
        {
            assert(wordCounts[sequenceInteger] > 0);

            if (wordReachedCutoff(sequenceInteger))
            {
                genomeIndex_t location =  basePairIndex - wordSize + 1;

                uint64_t sequenceIntegerLeft = 0;
                for (genomeIndex_t i = location - wordSize; i<location; i++) sequenceIntegerLeft = (sequenceIntegerLeft<<2) | gs.getInteger(i);
                wordHashLeft.addGenomeLocation(sequenceInteger, sequenceIntegerLeft, location);

                uint64_t sequenceIntegerRight = 0;
                for (genomeIndex_t i = location + wordSize; i<location + 2*wordSize; i++) sequenceIntegerRight = (sequenceIntegerRight<<2) | gs.getInteger(i);
                wordHashRight.addGenomeLocation(sequenceInteger, sequenceIntegerRight, location);

            }
            // #define WRITEHIGHREPEATLOCATIONS
#if defined(WRITEHIGHREPEATLOCATIONS)
            // for development purposes, this code prints to stdout the index word,
            // the location of this word, and the 30 bases at this location for
            // every case where the word occurs more than the cutoff number of
            // times - a high repeat pattern.
            if (wordReachedCutoff(sequenceInteger))
            {
                genomeIndex_t location =  basePairIndex - wordSize + 1;
                std::cerr << sequenceInteger
                          << " "
                          << location
                          << " ";
                for (genomeIndex_t i=0; i<2*wordSize; i++)
                {
                    std::cerr << gs[location + i];
                }
                std::cerr << "\n";
            }
#endif
            //
            // If this word has no entries, no more work to do.  This
            // will happen when either there were no entries, or when
            // the occurrence threshhold was exceeded.
            //
            if (hashindices[sequenceInteger] == hashindices[sequenceInteger+1])
                continue;

            // basePairIndex is pointing at the last character of this sequence, which is
            // why we adjust the position with '+ 1'
            wordpositions[hashindices[sequenceInteger+1] - wordCounts[sequenceInteger]] =
                basePairIndex - wordSize + 1;

            wordCounts[sequenceInteger]--;

            if (basePairIndex >= nextprogress)
            {
                int megabase = basePairIndex / 1000000;
                std::cout << megabase << " Mb\r" << std::flush;
                nextprogress += 10000000;
            }
        }
    } // end of for unsigned int i

#if defined(WRITEHIGHREPEATLOCATIONS)
    std::cerr << std::endl; //
#endif

    // sanity_check.storing_positions.txt
    std::cout << "DONE: Stored positions for non-highly repeated words" << std::endl << std::flush;

    delete [] wordCounts;
    wordCounts = NULL;
    //
    // To finish up, we write the malloc() version of wordpositions out to
    // disk.
    // We leave it hanging around in case the caller wants to use the word
    // index immediately after creating it.
    //
    int fd = ::open(wordPositionsFilename,O_CREAT|O_RDWR, 0666);
    if (fd<0)
    {
        perror("open:");
        exit(1);
    }
    //
    // write the header - magic cookie plus the version
    //
    ssize_t len, bytesWritten, wordPositionsLength;
    wordPositionsLength = sizeof(unsigned int) * (numberWordPositions);

    wordPositionMmapHeader *wpHeader = new wordPositionMmapHeader;
    wpHeader->cookie = UMWI_WP_COOKIE;
    wpHeader->version = UMWI_VERSION;
    wpHeader->numberWordPositions = numberWordPositions;
    wpHeader->wordPositionsOffset = sizeof(*wpHeader);
    wpHeader->wordReachedCutoffBitvectorOffset = wpHeader->wordPositionsOffset +
        wordPositionsLength;


    bytesWritten = write(fd, (char *) wpHeader, sizeof(*wpHeader));
    if (bytesWritten!=sizeof(*wpHeader))
    {
        std::cout << "FAILED to word position file header." << std::endl;
        exit(1);
    }

    bytesWritten = 0;
    // XXX more hacks... good grief...
    while (bytesWritten < wordPositionsLength)
    {
        /// length for ::write is size_t - 32bits?
        len = ::write(fd, ((char *) wordpositions) + bytesWritten, wordPositionsLength - bytesWritten);
        if (len==0)
        {
            perror("write:");
            exit(1);
        }
        bytesWritten += len;
    }
    len = write(fd, wordReachedCutoffBitvector, wordReachedCutoffBitvectorLength);
    if (len != (ssize_t) wordReachedCutoffBitvectorLength)
    {
        std::cout << "FAILED to write cutoff bit vector." << std::endl;
        exit(1);
    }
    ::close(fd);

    std::cout << "DONE: writing word positions arrays" << std::endl;
    fflush(stdout);
}


//
// for debugging:
//
void WordIndex::WriteOutWordPositions()
{
#if 0
    // XXX need to rewrite
    // if potential clashes: cannot guarantee hashindices[index] == index
    // therefore, cannot invoke this function since it will
    // replace hashindices[index] (i.e., sequenceInteger or wordAsInteger) with
    // starting_slot_for_(word)_positions

    for (unsigned int i = 0; i < wordsCountHashSize; i ++)
    {
        if (hashindices[i] == hashindices[i+1])
            continue;

        for (unsigned int j = hashindices[i]; j < hashindices[i+1]; j++)
            printf("%s %u\n", (const char *) IntegerToSeq(i), wordpositions[j]);
    } // end of i
#endif
}

//
// Given an existing word index, open it and map the hashindices and
// wordpositions arrays.
//
// Provided that there is enough RAM to cache them, there should be no
// disk reads after the first time this program has been run.
//
int WordIndex::open(GenomeSequence &reference)
{
    if (hashIndicesMemoryMap.open(hashIndicesFilename))
    {
        std::cerr << "Failed to open file '" << hashIndicesFilename << "'." << std::endl;
        return 1;
    }

    if (*((unsigned int *)hashIndicesMemoryMap.data + 0)!=UMWI_HI_COOKIE)
    {
        std::cerr << "file '" << hashIndicesFilename << "' has the wrong cookie number.  Expected "
                  << UMWI_HI_COOKIE << ", got " << *((unsigned int *)hashIndicesMemoryMap.data + 0) << "." << std::endl;
        return 1;
    }

    if (*((unsigned int *)hashIndicesMemoryMap.data + 1)!=UMWI_VERSION)
    {
        std::cerr << "file '" << hashIndicesFilename << "' has the wrong version number.  Expected "
                  << UMWI_VERSION << ", got " << *((unsigned int *)hashIndicesMemoryMap.data + 1) << "." << std::endl;
        return 1;
    }

    setWordSize(*((unsigned int *)hashIndicesMemoryMap.data + 2));

    //  assert(*((unsigned int *)hashIndicesMemoryMap.data + 3)==wordsCountHashSize);

    hashindices = (unsigned int *)hashIndicesMemoryMap.data + 3;

    //
    // map the word positions file
    //  - check cookie and version
    //  - set word position count
    //  - set word position array
    //  - set word reaches cutoff bitvector array
    //
    if (wordPositionsMemoryMap.open(wordPositionsFilename))
    {
        std::cerr << "Failed to open file '" << wordPositionsFilename << "'." << std::endl;
        return 1;
    }

    wordPositionMmapHeader *wpHeader = (wordPositionMmapHeader *) wordPositionsMemoryMap.data;

    if (wpHeader->cookie != UMWI_WP_COOKIE)
    {
        std::cerr << "file '" << wordPositionsFilename << "' has the wrong cookie.  Expected "
                  << UMWI_WP_COOKIE << ", got " << wpHeader->cookie << "." << std::endl;
        return 1;
    }

    if (wpHeader->version != UMWI_VERSION)
    {
        std::cerr << "file '" << wordPositionsFilename << "' has the wrong version number.  Expected "
                  << UMWI_VERSION << ", got " << wpHeader->version << "." << std::endl;
        return 1;
    }

    numberWordPositions = wpHeader->numberWordPositions;

    wordpositions = (unsigned int *)((char *) wordPositionsMemoryMap.data + wpHeader->wordPositionsOffset);

    wordReachedCutoffBitvector = (unsigned char *)((char *) wordPositionsMemoryMap.data + wpHeader->wordReachedCutoffBitvectorOffset);

#if 0
    testHashIndices(reference);
    exit(1);
#endif

    return 0;
}

//
// Top level method for creating the word index for the given
// reference genome.
//
// To preserve mapping sensitivity, this method also populates
// the left and right high repeat hashes.  The right hash is formed
// from the key {highRepeatWord, nextWordInRead} and the left hash
// is formed from the key {previosWordInRead, highRepeatWord}.
//
//
int WordIndex::create(GenomeSequence &reference, WordHash &wordHashLeft, WordHash &wordHashRight, int wordSize)
{
    Timing totalCreateTime;
    Timing createTime;

    setWordSize(wordSize);  // set word size, as well as hash table size and mask

    totalCreateTime.start();

    createTime.start();

    // populate the hashindices array
    WordsCountHash(reference);

    createTime.end();

    std::cout << "Elapsed time to create hashindices: " << createTime.interval()/60 << " minutes." << std::endl;

    createTime.start();

    // Populate the word count array
    IndexPositionsForLowRepeatWords();

    createTime.end();

    std::cout << "Elapsed time to populate the wordcount array: " << createTime.interval()/60 << " minutes." << std::endl;

    createTime.start();

    // populate the word positions array
    StoreWordPositions(reference, wordHashLeft, wordHashRight);

    createTime.end();

    std::cout << "Elapsed time to populate word positions array: " << createTime.interval()/60 << " minutes." << std::endl;

    totalCreateTime.end();

    std::cout << "Total create time: " << totalCreateTime.interval()/60 << " minutes.\n" << std::endl;

    return 0;
}

//
// populate the various filenames we're using.
//
void WordIndex::setFilenames(const String &base)
{
    baseFilename = base;

    wordPositionsFilename = baseFilename + ".umwiwp";
    hashIndicesFilename = baseFilename + ".umwihi";
}

//
// Given a word, find the total number of possible matches for
// this word, including the word itself plus all single base
// mutations.
//
int WordIndex::getMatchCountWithMutations(wordInteger_t word)
{
    if (word==INVALID_WORDINDEX) return 0;

    int count = getMatchCount(word);    // start with no mutations

    wordInteger_t mask1 = 1;
    wordInteger_t mask2 = 2;
    wordInteger_t mask3 = 3;
    for (unsigned int i=0; i<wordSize; mask1<<=2, mask2<<=2, mask3<<=2, i++)
    {
        count += getMatchCount(word ^ mask1);
        count += getMatchCount(word ^ mask2);
        count += getMatchCount(word ^ mask3);
    }
    return count;
}

//
// Given a single word, check all locations of that word and ensure
// that the matches correctly point at the indexed word.
//
int WordIndex::testHashIndex(GenomeSequence &gs, wordInteger_t word)
{
    int rc = 0;
    unsigned int wpIndex;

    if (hashindices[word] ==  hashindices[word+1]) return 0; // nothing to do

    assert(!wordReachedCutoff(word));

    String bases = IntegerToBPSeq(word);

    for (wpIndex = hashindices[word]; wpIndex < hashindices[word+1]; wpIndex++)
    {
        int i;
        unsigned int genomeIndex = wordpositions[wpIndex];
        for (i = 0; i<bases.Length(); i++)
        {
            if (gs[genomeIndex + i] != bases[i])
            {
                rc = 1;
                break;
            }
        }
        if (i<bases.Length())
        {
            std::cerr << "wpIndex[" << wpIndex - hashindices[word]
                      << "] = " << wpIndex << ", word 0x" << std::hex << word << std::dec << ": expect: " << bases << " got ";
            for (i = 0; i<bases.Length(); i++)
            {
                if (gs[genomeIndex + i] == bases[i])
                {
                    std::cerr << gs[genomeIndex + i];
                }
                else
                {
                    std::cerr << tolower(gs[genomeIndex + i]);
                }
            }
            //          gs.printNearbyWords(genomeIndex, 60, bases);
            std::cerr << std::endl;
        }
    }

    return rc;
}

//
// testHashIndices loops over all possible words, and verifies that
// the matches it has are actually in the genome sequence.  This
// should validate the hashindices and wordpositions arrays.
//
void WordIndex::testHashIndices(GenomeSequence &gs)
{

    int rc = 0;
    std::cout << "checking all possible words:" << std::endl;
    for (uint64_t word=0; word<wordsCountHashSize; word++)
    {
        rc += testHashIndex(gs, word);
    }

    std::cout << rc << " words were bad." << std::endl;

    std::cout << "now testing for words listed in more than one hash index." << std::endl;

    std::map<genomeIndex_t, unsigned int>    index2word;

    int wtf = 0;
    // XXX long long unsigned int shuts up compiler warnings... should use uint64_t
    for (long long unsigned int word=0; word<wordsCountHashSize; word++)
    {
        for (unsigned int wpIndex = hashindices[word]; wpIndex < hashindices[word+1]; wpIndex++)
        {
            unsigned int genomeIndex = wordpositions[wpIndex];
            if (index2word.find(genomeIndex) == index2word.end())
            {
                index2word[genomeIndex] = word;
            }
            else
            {
                if (index2word[genomeIndex] == word)
                {
                    // wound up mapping a location twice in the same word - bizarre
                    std::cout << "genome position " << genomeIndex << " duplicated for word " << word << std::endl;
                }
                else
                {
                    // wound up mapping a location once for one word, and again for another
                    std::cout << "genome position " << genomeIndex << " duplicated for word "
                              << word << " and " << index2word[genomeIndex] << "!" << std::endl;
                    wtf++;
                }
            }
        }
    }
    std::cout << "impossible condition occurs " << wtf << " times." << std::endl;
}

void WordIndex::whatWordPointsAt(genomeIndex_t g)
{
    for (uint64_t word=0; word<wordsCountHashSize; word++)
    {
        for (unsigned int wpIndex = hashindices[word]; wpIndex < hashindices[word+1]; wpIndex++)
        {
            unsigned int genomeIndex = wordpositions[wpIndex];
            if (g == genomeIndex)
            {
                std::cout << "word " << std::hex <<
                    word <<
                    std::dec <<
                    " maps to genome location " << g << std::endl;
            }
        }
    }
}

//
// XXX unfinished - need to merge wordSize * 4 different lists
// of positions together into the positions vector.
//
bool WordIndex::getWordPositions(wordInteger_t word, std::vector<genomeIndex_t> &positions)
{
    int count = 0;
    genomeIndex_t *candidates = NULL;

    wordInteger_t mask1 = 1;
    wordInteger_t mask2 = 2;
    wordInteger_t mask3 = 3;

    positions.clear();

    // each base is represented as two bits in the index word.
    // there are 4 combinations for each base.
    // the total number of sets of positions is therefore 4 x wordSize.
    std::vector<int> counts(4 * wordSize);
    std::vector<int> offsets(4 * wordSize);

    int offset = 0;
    for (unsigned int i=0; i<wordSize; mask1<<=2, mask2<<=2, mask3<<=2, i++)
    {

        //
        // NB: std::copy might be better than vector<T>::insert
        //
        // For each word we are interested in, push back a position, a count,
        // and a copy of the word positions.  When done, we will do a recursive
        // inplace_merge, which will tend to minimize data movement.
        //
#define PUSH_BACK_POSITIONS(mask)                                       \
        count = hashindices[(word ^ mask) + 1] - hashindices[(word ^ mask)]; \
        if(count) {                                                     \
            candidates = &wordpositions[hashindices[(word ^ mask)]];    \
            positions.insert(positions.end(), candidates, candidates+count); \
            counts.push_back(count);                                    \
            offsets.push_back(offset);                                  \
            offset += count;                                            \
        }

        PUSH_BACK_POSITIONS(0);
        PUSH_BACK_POSITIONS(mask1);
        PUSH_BACK_POSITIONS(mask2);
        PUSH_BACK_POSITIONS(mask3);

    }
    //
    // positions.size() may be 0 here (no matching positions for this index word)
    // counts and offsets are the index into the vector positions of all the
    // sorted ranges.
    // Now use inplace_merge to put the vector positions into ascending order.
    //
    inplace_merge(offsets, counts, 0, offsets.size(), positions);
    return false;
}

