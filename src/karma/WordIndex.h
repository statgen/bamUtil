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

#ifndef __WORD_INDEX_H__
#define __WORD_INDEX_H__

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>

#include "GenomeSequence.h"
#include "IntHash.h"
#include "MemoryMap.h"
#include "Sort.h"
#include "StringArray.h"
#include "WordHash.h"


struct wordPositionMmapHeader
{
    uint32_t    cookie;
    uint32_t    version;
    uint32_t    numberWordPositions;
    size_t wordPositionsOffset;    // from beginning of file
    size_t wordReachedCutoffBitvectorOffset;
};


#if !defined(HAVE_WORDINTEGER_T)
#define HAVE_WORDINTEGER_T
typedef uint32_t wordInteger_t;
typedef uint32_t wordIndex_t;
#define INVALID_WORDINDEX UINT32_MAX
#endif

class WordIndex
{
 protected:
    int debugFlag;
    int verboseFlag;

    static const unsigned int UMWI_VERSION = 20090617U;  // YYYYMMDD of last change to file layout
    String      baseFilename;
    String IntegerToBPSeq(wordInteger_t word);

 public:
    WordIndex();
    ~WordIndex();

    void setVerbose(bool v)
    {
        verboseFlag = v;
    }
    void setDebug(bool d)
    {
        debugFlag = d;
    }

    // UMWI -> UM Word Index
    static const unsigned int UMWI_HI_COOKIE = 0xe227a738;  // some unique random cookie ID
    static const unsigned int UMWI_WP_COOKIE = 0xd242d62c;  // some unique random cookie ID


    void        setFilenames(const String &base);
    String      hashIndicesFilename;
    String      wordPositionsFilename;
    MemoryMap   hashIndicesMemoryMap;
    MemoryMap   wordPositionsMemoryMap;

    unsigned int occurrenceCutoff;       // maximum mismatch count (e.g. 4)

    uint32_t wordSize;
    uint64_t wordsCountHashSize;
    uint64_t wordsCountIndexMask;

    unsigned int numberWords; // including highly repeated words;

    unsigned int numberWordPositions; // for non-highly repeated words;

    void setWordSize(int wordSize)
    {
        this->wordSize = wordSize;
        wordsCountHashSize = (1ULL << (wordSize << 1));
        wordsCountIndexMask = wordsCountHashSize - 1;
        wordReachedCutoffBitvectorLength = (wordsCountHashSize + 7) / 8;
    }

#define WORDCOUNT_MAX UINT32_MAX
    // temporarily stores sequenceInteger (or wordAsInteger) after WordsCountHash();
    // changes to starting_slot_for_(word)_positions.  Used during only ::create().
    unsigned int *  wordCounts;

    unsigned char *wordReachedCutoffBitvector;
    unsigned int wordReachedCutoffBitvectorLength;  // size in bytes, not bits

 public:
    unsigned int * hashindices; // given a word index, point to a list of words in wordpositions

    genomeIndex_t   *wordpositions; // an array of lists of words

    int open(GenomeSequence &reference);
    void prefetch()
    {
        wordPositionsMemoryMap.prefetch();
        hashIndicesMemoryMap.prefetch();
    }

    // use reference to create a short word index (15 bases) and also
    // a longer 30 base hash for breaking high repeat cases.
    int create(GenomeSequence &reference, WordHash &wordHashLeft, WordHash &wordHashRight, int wordSize = 15);

    void setOccurrenceCutoff(int co)
    {
        occurrenceCutoff = co;
    }

    //
    // true->the given word appears in the genome too many
    // times, so it isn't included in the word positions array.
    //
    bool wordReachedCutoff(unsigned int word);
    //
    // Fast check to see how many repeats there are of the given word
    //
    int getMatchCount(wordInteger_t word);
    //
    // slower check to see how many repeats there are of the given
    // word, but including all one character mutations.
    //
    int getMatchCountWithMutations(wordInteger_t word);

    int testHashIndex(GenomeSequence &gs, wordInteger_t word);
    void testHashIndices(GenomeSequence &gs);

    // given a word and a target search area, return a match within there
    unsigned int   findWordPosition(wordInteger_t word, genomeIndex_t filter, uint32_t width);

    inline genomeIndex_t   getMatchPosition(
                                            unsigned int wordPositionIndex,
                                            int wordStartPosition
                                            ) const
    {
        return  wordpositions[wordPositionIndex] - wordStartPosition;
    }

    const String &getHashIndicesFilename() const
    {
        return hashIndicesFilename;
    }


 protected:
    bool AllocateWordsCountHashMemory();

    // only hashes the first instance of each wordSize sequence;
    void WordsCountHash(GenomeSequence & gs);
    void WriteWordsCountHash(const String & prefix);

    void IndexPositionsForLowRepeatWords();

    // XXX not implemented:
    int GetOccurrenceCutoff(double pctWordsToKeep);

    void IndexPositionsForLowRepeatWords(double pctWordsToKeep);

    // wordCounts and wordRunningCounts will be reset to zero in this function;
    void StoreWordPositions(GenomeSequence & gs, WordHash &wordHashLeft, WordHash &wordHashRight);

    // for debug:
    void WriteOutWordPositions();
    void whatWordPointsAt(genomeIndex_t);

    bool getWordPositions(wordInteger_t word, std::vector<genomeIndex_t> &positions);
};

std::ostream &operator << (std::ostream &stream, WordIndex &);

inline bool WordIndex::wordReachedCutoff(wordInteger_t word)
{
    assert(word>>3 < wordReachedCutoffBitvectorLength);
    return (wordReachedCutoffBitvector[word>>3] >> (word&0x7)) & 0x1;
}

//
// Note: does not check for INVALID_WORDINDEX
//
inline int WordIndex::getMatchCount(wordInteger_t word)
{
    return hashindices[word+1] - hashindices[word];
}


//
// for a given word, find if there is a match within window in the range
//  filter - width .. filter + width.
//
inline unsigned int WordIndex::findWordPosition(wordInteger_t word, genomeIndex_t filter, uint32_t width)
{
    wordIndex_t lo = hashindices[word];      // inclusive lower bound
    wordIndex_t hi = hashindices[word+ 1];   // exclusive
    wordIndex_t mid;

    // the binary search below works for most values, but not
    // when both lo and hi are 0, so may as well leave here:
    if (lo==hi) return INVALID_WORDINDEX;

    hi--;   // now an inclusive upper bound
    while (lo <= hi)
    {
        mid = lo + (hi - lo) / 2;

        if (wordpositions[mid] > filter + width)
        {
            // the following test is because the lo<=hi test only works
            // for signed index values, but we're unsigned, so hi, rather
            // than being -1, is MAX_INT, so the loop continues forever.
            if (mid==0) return INVALID_WORDINDEX;
            hi = mid - 1;
        }
        else if (wordpositions[mid] < filter - width)
        {
            lo = mid + 1;
        }
        else
        {
            return mid;
        }
    }
    return INVALID_WORDINDEX;
}

#endif
