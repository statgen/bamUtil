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

#include "ReadIndexer.h"
#include "GreedyTupleAligner.h"
#include "TrimSequence.h"
#include "SmithWaterman.h"
#include <algorithm>    // for sort
#include <iostream>

int ReadIndexer::Word2Integer(const char *read, int offset, int &wordNLocation)
{
    //    assert(word.size() + index >= wordIndex->wordSize && index >= 0);

    int wordinteger = 0;
    wordNLocation = -1;

    for (unsigned int i = offset; i < offset + wordIndex->wordSize; i++)
    {
        wordinteger <<= 2;
        int crtbaseinteger = GenomeSequence::base2int[(int) read[i]];

        switch (crtbaseinteger)
        {
        case GenomeSequence::baseXIndex:
            // policy: if a single index word in a read contains
            // a single unknown character, we currently consider
            // it different from an 'N', and discard the index.
            // The read can still be mapped, just not on the
            // basis of this particular index word position.
            return INVALID_WORDINDEX;
        case GenomeSequence::baseNIndex:
            // policy: indexing allows for a single N character to exist
            // in each word.  If more than one exists, we can't readily
            // code to mutate all possible combinations, so we'll bail in
            // that case.
            if (wordNLocation<0)
            {
                wordNLocation = i - offset;
                crtbaseinteger = 0;         // knowing this is 0, we'll visit all edits of this base later
            }
            else
            {
                return INVALID_WORDINDEX;
            }
            // drop through to default handling
        default:
            wordinteger |= crtbaseinteger;
        }
    }

    return wordinteger;
}

ReadIndexer::ReadIndexer(MapperUserOption &m) : mapperOptions(m)
{
    readLength = -1;
    isForward = true;
    gs = NULL;
    wordIndex = NULL;
    checkEdits = false;
    useGapped = false;
}

struct wordsCompareClass {
    bool operator() (const ReadIndexer::Word& w1, const ReadIndexer::Word& w2)
    {
        return (w1.sortKey < w2.sortKey);
    }
} wordsCompareInstance;

//
// Set overall index matching strategy for this read.
//
// It is different for forward versus backward, so we
// have to do it in this class.
//
// The end result is to populate wordInets, wordNLocations and wordPositions
// with the words that we will be using.  In addition, it sets the matchStrategy.
//
// XXX
//  what are the desired end results?
//  how and when will we check the secondary hash?
//  how many basic iteration cases are there? edit/noedit/what?
//
//  Geez...
//
//  this needs to be driven by how we choose to assemble and
//  iterate over the candidates.
//
// @return bool true: failed to set index strategy; false: succeed
bool ReadIndexer::setIndexStrategy()
{
    checkedPositions.Clear();       // clear hash of genome positions we have already checked.

    // the loops below are not zero size friendly...
    if (read.size() < wordIndex->wordSize) return true;

    words.resize(0);

    //
    // by default, mapperOptions.checkIndexWordMutations is false, but
    // it can be overridden for testing.  Typically this would be a
    // test to see if we can increase sensitivity at the expense of
    // additional runtime.
    //
    checkEdits = mapperOptions.checkIndexWordMutations;

    // collect all the possible words into a vector.
    //
    // later, we may wish to collect only "enough" to ensure
    // that the rest of this method yields reasonable results
    // for most reads.
    //
    int lowRepeatWords = 0;
    for (unsigned int offset = 0; offset < read.size() - wordIndex->wordSize; offset += wordIndex->wordSize)
    {
        struct Word word;

        if (isForward) word.position = offset;
        else word.position =
                 readLength - offset - wordIndex->wordSize;

        word.wordInt =  Word2Integer(read.c_str(), word.position, word.NLocation);

        //
        // words with more than one 'N' in them are difficult
        // to deal with, so we drop them here.
        //
        if (word.wordInt==INVALID_WORDINDEX) continue;

        // force high repeat words to the back of the list
        word.sortKey = UINT32_MAX;
        if (!wordIndex->wordReachedCutoff(word.wordInt))
        {
            word.sortKey = wordIndex->getMatchCount(word.wordInt);
            lowRepeatWords ++;
        }
        words.push_back(word);
    }

    wordInts.resize(0);
    wordNLocations.resize(0);
    wordPositions.resize(0);

    //
    // right now, this is slow, but make it work first.
    //
    // This sorts the words list by the number of genome
    // match positions, low to high.  The goal is to cherry
    // pick and take the most unique indices.
    //
    std::sort(words.begin(), words.end(), wordsCompareInstance);
    //
    // What we ideally want to do is populate wordPositions with the
    // smallest number of non-repeat words possible.
    //
    // Here we cherry pick - if we get 3 non-repeat index words,
    // we'll call it good enough.
    //
    // XXX paramaterize magic number "3" - I think 3 is fine,
    // Goncalo prefers 4.  I think what I want is 3 by default,
    // but a parameter to adjust to 4 so we can test performance
    // and sensitivity on real data.
    //
    // NB: An index word with a single 'N' in it is actually only
    // a little bit worse than a clean index word, so for purposes
    // of decision making in this method, we'll treat it as being
    // the same as a clean index word.  Of course, the matching code
    // will need to iterate over the 4 base values.
    //
#define MIN_INDEX_COUNT 2

    checkEdits = false;
    unsigned int whichWord = 0;
    while (wordInts.size() < mapperOptions.maximumIndexCount &&   // while not enough words
           whichWord < words.size()              // while more exist
           )
    {
        wordInts.push_back(words[whichWord].wordInt);
        wordNLocations.push_back(words[whichWord].NLocation);
        wordPositions.push_back(words[whichWord].position);
        whichWord++;
    }

    if (wordInts.size()==mapperOptions.maximumIndexCount)
    {
        // hey hey hey, we have enough!
        return false;
    }

    if (wordInts.size()>=mapperOptions.minimumIndexCount)
    {
        // we have enough clean index words, but we need
        // to check all edits on them.  Before we finish,
        // ensure we have the minimum number of index
        // words:
        if (wordInts.size() > mapperOptions.minimumIndexCount)
        {
            wordInts.resize(mapperOptions.minimumIndexCount);
            wordNLocations.resize(mapperOptions.minimumIndexCount);
        }
        checkEdits = true;
        return false;
    }

    //
    // If we got here, it means that we have too few
    // useful index words to get a map.  This happens
    // either when a read is too short, or there are too many
    // "bad" index words (each index word may contain up to
    // 1 invalid base - typically an 'N' for base space reads).
    //

    return true;

}

//
// for color space or base space, set up the read and quality
// values.
//
// color space has different handling due to its encoding being
// symmetric in part.
//
// for base space only, first need to do quality trimming (-q # option).
//
//
void ReadIndexer::setReadAndQuality(const char *readArg, int readLengthArg, const char *qualityArg)
{
    this->readLength = readLengthArg;

    // default is the whole read is good:
    int firstGoodBase = 0;          // closed range
    int lastGoodBase = readLength;  // open range

    // a lot to bookkeep, but we have no choice:
    leftTrimmedBases.clear();
    rightTrimmedBases.clear();
    leftTrimmedQualities.clear();
    rightTrimmedQualities.clear();

    // but sometimes we trim on quality
    if (!isdigit(*readArg) && mapperOptions.qualityTrim)
    {
        // For base space, when user asks for it, we trim based on quality.
        //
        // When we trim, we have to keep track of the trimmed bits so we can
        // print them as soft clipped portions in the read.
        //
        binaryQuality.clear();
        for (int i=0; i<readLengthArg; i++) binaryQuality.push_back(qualityArg[i] - '!');
        firstGoodBase = trimSequence(binaryQuality, mapperOptions.qualityTrim, true) - binaryQuality.begin();
        lastGoodBase = trimSequence(binaryQuality, mapperOptions.qualityTrim, false) - binaryQuality.begin() + 1;

        // after trimming, check for sane results.
        // One case is there are no good bases anywhere (or no run of
        // more than 4 bases with good quality).
        // Another case is when there are between 1 and howManyValues
        // good bases followed by bad bases... we'll wind up on the
        // right hand side of the very short good region.
        if (lastGoodBase <= firstGoodBase)
        {
            leftTrimmedBases = readArg;
            leftTrimmedQualities = qualityArg;
            this->readLength = 0;
            this->read = "";
            this->phredQuality = "";
            return;
        }

        this->readLength = lastGoodBase - firstGoodBase - 1;    // lastGoodBase actually points one after the last

        leftTrimmedBases.assign(readArg, firstGoodBase);
        rightTrimmedBases.assign(readArg + lastGoodBase - 1, readLengthArg - lastGoodBase + 1);

        leftTrimmedQualities.assign(qualityArg, firstGoodBase);
        rightTrimmedQualities.assign(qualityArg + lastGoodBase - 1, readLengthArg - lastGoodBase + 1);

        // sanity check
        assert(leftTrimmedBases.size() + rightTrimmedBases.size() + this->readLength == (unsigned) readLengthArg);

        //
        // Now reverse the trimmed fragments according to technology
        //
        // We could mix this in below, but I'd rather keep it clean and clear as
        // practical here.
        //

        if (!isForward)
        {
            // qualities get swapped regardless of technology
            reverse(leftTrimmedQualities.begin(), leftTrimmedQualities.end());
            reverse(rightTrimmedQualities.begin(), rightTrimmedQualities.end());
            // bases/colors are swapped differently
            if (isdigit(read[0]))
            {
                // Color Space:
                std::reverse(leftTrimmedBases.begin(), leftTrimmedBases.end());
                std::reverse(rightTrimmedBases.begin(), rightTrimmedBases.end());
            }
            else
            {
                // Base space:
                gs->getReverseRead(leftTrimmedBases);
                gs->getReverseRead(rightTrimmedBases);
            }
        }
    }

    // good bases and quality
    std::string goodBases(readArg + firstGoodBase, this->readLength);
    std::string goodQuality(qualityArg + firstGoodBase, this->readLength);

    //
    // now all indexing should be done in the range of [firstGoodBase, lastGoodBase].
    // all reference to base or quality values should be to goodBases and goodQuality,
    // respectively.
    //

    if (isForward)
    {
        this->read = goodBases;
        phredQuality = goodQuality;

        packedRead.set(this->read.c_str(), 0);
        packedReadWithPad.set(this->read.c_str(), 1);
        binaryQuality.clear();
        for (size_t i=0; i<goodQuality.size(); i++) binaryQuality.push_back(goodQuality[i] - '!');
    }
    else
    {

        //
        // Handle reverse case.
        //
        // according to SOLiD_Dibase_Sequencing_and_Color_Space_Analysis.pdf, forward and
        // reverse complement colors do not change, as a consequence, if we detect
        // an ABI SOLiD read here, don't complement the colors...
        //
        // NB: here the first symbol is the symbol following the previously truncated
        // primer base.  For ABI SOLiD, it should always be a digit.
        //
        // NB: although ABI SOLiD doesn't have color values swapped, QUALITY VALUES
        // MUST STILL BE SWAPPED.
        //
        this->read.clear();
        binaryQuality.clear();
        if (isdigit(readArg[0]))
        {
            // ABI SOLiD
            for (int i = this->readLength - 1; i >= 0; i --)
            {
                // don't complement the color:
                this->read.push_back(goodBases[i]);
                binaryQuality.push_back(goodQuality[i] - '!');
            }
        }
        else
        {
            // base space
            for (int i = this->readLength - 1; i >= 0; i --)
            {
                // complement the base:
                this->read.push_back(GenomeSequence::base2complement[(int) goodBases[i]]);
                binaryQuality.push_back(goodQuality[i] - '!');
            }
        }
        packedRead.set(this->read.c_str(), 0);
        packedReadWithPad.set(this->read.c_str(), 1);

        // also reverse phred encoded quality - we already reversed
        // the binary qualities above
        phredQuality = goodQuality;
        reverse(phredQuality.begin(), phredQuality.end());
    }
}

void ReadIndexer::dump()
{
    std::cout << (isForward ? "Forward" : "Backward") << " strand index dump:\n";
    std::cout << "   Read: " << read << "\n";
    std::cout << "Quality: " << phredQuality << "\n";
    std::cout << "         Read length : " << read.size() << "\n";
    std::cout << "     index word count: " << wordInts.size() << "\n";
    for (uint32_t whichWord = 0; whichWord < wordInts.size(); whichWord++)
    {
        std::cout << "Word " << whichWord << ": " << std::hex << wordInts[whichWord] << std::dec
                  << " starts at " << wordPositions[whichWord];
        if (wordInts[whichWord]!= INVALID_WORDINDEX && wordNLocations[whichWord]>-1)
        {
            std::cout << " and has exactly one N at offset " << wordNLocations[whichWord];
        }
        std::cout << std::endl;
    }
}


//
// getSumQSW gets the quality and cigarString of a read at a given
// genome position using the Smith Waterman algorithm.
//
// This is expected to be slow code.  We do a banded approach so that
// we don't have to clear or populate the full H matrix.  There are
// 'flat' implementations that save storage by normalizing the band
// so that it is horizontal.
//
// The overall strategy is roughly this:  since we know that the bases
// in the index word have no indels, we don't need to run smith waterman
// there.  We do need to get the sumQ, of course.  To save time and
// space, we run Smith Waterman twice - once for the bases from just after
// the index word to the end of the read, and again for the bases just
// prior to the index word, through the first base of the read.
//
// This brings up the issue of running Smith Waterman "forwards" versus
// "backwards" - since the bases after the index word are in the forward
// direction, things run smoothly.  However, the bases just before the
// index word differ in the sense that if you run Smith Waterman from the
// beginning of the read through the base before the index word, you lose
// the ability to handle the same number of indels as you do in the string
// of bases following the index word.  To keep this balanced, the
// Smith Waterman algorithm used here allows for being run in reverse, which
// simply means running from the end of the string to the head of the string.
//
// There are some nuances to this that are difficult to deal with
// even in the test code in libcsg/SmithWaterman.cpp.
//
// One code path that is messy is when the read has 1 or more deleted
// bases in it.  In those cases, you don't simply want to run the
// SW algorithm against the "read length" number of bases, because if
// you do, there will be bases at the end of the read that can't be
// mapped to the reference, if you specify the reference to be the same
// length as the read.  Therefore, as a rule, you need to expand the
// reference length by the number of bases that you want to allow to be
// deleted in a read.  This number is mapperOptions.smithWatermanBandSize
// minus one.  I left out the -1 below to give myself wiggle room.  Who
// knows if that is worthwhile...
//
// Most of the other code paths are dealt with reasonably well in the test
// code.
//
//
int ReadIndexer::getSumQSW(
                           genomeIndex_t matchPosition,
                           int &mismatchCount,
                           int bestMatchSumQ,
                           int whichWord)
{
    mismatchCount = 0;
    // all this work is to return sumQ:
    int sumQ = 0;

    //
    // be careful with uninitialized data.  The constructor clears the entire
    // matrix H, which is slow.  Do this only once each time setAllowedInsertDelete
    // is changed.
    //
    static SmithWaterman<1024, 1024, uint16_t, std::string , GenomeSequence , std::string, uint32_t, genomeIndex_t > SW;

    SW.clear(); // only resets member variables, not the matrix H - use clearH() for that.

    //
    // First, we wish to do a Smith Waterman alignment from the
    // base following the index word through the end of the string.
    //
    SW.setRead(&read);
    SW.setReadQuality(&phredQuality);
    SW.setReference(gs);
    SW.setAllowedInsertDelete(mapperOptions.smithWatermanBandSize);

    int readAndReferenceLength = read.size() - (wordPositions[whichWord] + wordIndex->wordSize);
    SW.setDirection(+1);        // SW is statically allocated, always need to set this
    if (readAndReferenceLength > 0)
    {


        // point to the base just after the index word:
        SW.setReadOffset(wordPositions[whichWord] + wordIndex->wordSize);
        SW.setReferenceOffset(matchPosition + wordPositions[whichWord] + wordIndex->wordSize);

        // set the match length to all bases following the index word:
        SW.setReadLength(readAndReferenceLength);

        // the addition of the band size allows a read to have up
        // to bandsize-1 deletes and still match the remainder bases
        // against the reference without calling them soft clips.
        SW.setReferenceLength(readAndReferenceLength + mapperOptions.smithWatermanBandSize);

        // SLOW: compute the MxN matrix (ok, not really - it is banded)
        SW.populateH();
#if 1
        // XXX WRONG:
        if (SW.getSoftClipCount() > ((int) read.size() > 40 ? (int) read.size()/20 : 0))
        {
            return INT_MAX;
        }
#endif
        // slightly faster: populate the alignment vector
        SW.populateAlignment();

        // slightly faster yet: calculate the sumQ
        sumQ += SW.getSumQ();

        if (mapperOptions.debug)
        {
            SW.printH();
        }
    }
    else
    {
        SW.clearAlignment();
    }

    //
    // check for mismatches inside the index word portion of the read.
    //
    for (genomeIndex_t i = wordPositions[whichWord]; i < wordPositions[whichWord] + wordIndex->wordSize; i++)
    {
        if (read[i] != (*gs)[matchPosition + i]) sumQ += phredQuality[i];
    }


    //
    // Now do SW backwards from the base immediately preceding the
    // index word to the start of the read:
    //

    readAndReferenceLength = wordPositions[whichWord];
    if (readAndReferenceLength > 0)
    {
        SW.maxCostPosition.first = SW.maxCostPosition.second = 0;

        SW.setDirection(-1);        // SW is statically allocated, always need to set this

        // point to the first base:
        SW.setReadOffset(0);
        SW.setReferenceOffset(matchPosition - mapperOptions.smithWatermanBandSize);

        // set the match length to all bases up to the index word:
        SW.setReadLength(readAndReferenceLength);

        // allow bandsize - 1 bases to be deleted in the read (see
        // related comments in similar code)
        SW.setReferenceLength(readAndReferenceLength + mapperOptions.smithWatermanBandSize);

        // SLOW: compute the MxN matrix (ok, not really - it is banded)
        SW.populateH();
#if 1
        //
        // this checks to see how much of the read we are willing to clip
        // off.  If the read is longer than 40 bases, allow up to 5% soft
        // clipping, otherwise allow no soft clipping.
        //
        if (SW.getSoftClipCount() > ((int) read.size() > 40 ? (int) read.size()/20 : 0))
        {
            return INT_MAX;
        }
#endif
        // slightly faster: populate the alignment vector
        SW.populateAlignment();

        // slightly faster yet: calculate the sumQ
        sumQ += SW.getSumQ();

        if (mapperOptions.debug)
        {
            SW.debugPrint();
        }
    }
    else
    {
        SW.clearAlignment();
    }

    return sumQ;
}


//
// getCigarString gets the cigarString of a read at a given
// genome position using the Smith Waterman algorithm.
//
// XXX rename cigarRollberBackward, since it is now an
// argument to this function (since we now return it by
// reference, so the best match print routines can correctly
// print matches/mismatches against the reference in the face
// of indels).
//
int ReadIndexer::getCigarRoller(
                                genomeIndex_t matchPosition,
                                int whichWord,
                                int &matchPositionOffset,
                                CigarRoller &cigarRollerBackward)
{
    CigarRoller cigarRollerForward;
    cigarRollerBackward.clear();

    static SmithWaterman<1024, 1024, uint16_t, std::string , GenomeSequence , std::string, uint32_t, genomeIndex_t > SW;

    SW.clear(); // only resets member variables, not the matrix H - use clearH() for that.

    //
    // First, we wish to do a Smith Waterman alignment from the
    // base following the index word through the end of the string.
    //
    SW.setRead(&read);
    SW.setReadQuality(&phredQuality);
    SW.setReference(gs);
    SW.setAllowedInsertDelete(mapperOptions.smithWatermanBandSize);

    int readAndReferenceLength = read.size() - (wordPositions[whichWord] + wordIndex->wordSize);
    SW.setDirection(+1);        // SW is statically allocated, always need to set this
    if (readAndReferenceLength > 0)
    {


        // point to the base just after the index word:
        SW.setReadOffset(wordPositions[whichWord] + wordIndex->wordSize);

        SW.setReferenceOffset(matchPosition + wordPositions[whichWord] + wordIndex->wordSize);

        // set the match length to all bases following the index word:
        SW.setReadLength(readAndReferenceLength);

        // here, the match length is expanded to allow for the max number
        // of deletes in the read, otherwise, the read CIGAR string will
        // show up as #S (some number of soft clips).  The smithWatermanBandSize
        // is (I think) 1 larger than the max number of indels that will
        // be detected.
        SW.setReferenceLength(readAndReferenceLength + mapperOptions.smithWatermanBandSize);

        // SLOW: compute the MxN matrix (ok, not really - it is banded)
        SW.populateH();
#if 1
        if (SW.getSoftClipCount() > ((int) read.size() > 40 ? (int) read.size()/20 : 0))
        {
            return INT_MAX;
        }
#endif
        // slightly faster: populate the alignment vector
        SW.populateAlignment();

        if (mapperOptions.debug)
        {
            SW.printH();
        }
    }
    else
    {
        SW.clearAlignment();
    }

    // This depends on direction being set to +1
    SW.rollCigar(cigarRollerForward);

    //
    // Now do SW backwards from the base immediately preceding the
    // index word to the start of the read:
    //

    readAndReferenceLength = wordPositions[whichWord];
    SW.setDirection(-1);        // SW is statically allocated, always need to set this
    if (readAndReferenceLength > 0)
    {
        SW.maxCostPosition.first = SW.maxCostPosition.second = 0;


        // point to the first base:
        SW.setReadOffset(0);

        // here we are going backwards from the end of the two
        // strings to the front, but we want to allow for indels in
        // the read to continue matching bases against the reference,
        // so we need to allow additional bases in the reference for that.
        SW.setReferenceOffset(matchPosition - mapperOptions.smithWatermanBandSize);

        // set the match length to all bases up to the index word:
        SW.setReadLength(readAndReferenceLength);
        SW.setReferenceLength(readAndReferenceLength + mapperOptions.smithWatermanBandSize);

        // SLOW: compute the MxN matrix (ok, not really - it is banded)
        SW.populateH();
#if 1
        //
        // this checks to see how much of the read we are willing to clip
        // off.  If the read is longer than 40 bases, allow up to 5% soft
        // clipping, otherwise allow no soft clipping.
        //
        // XXX FIX THIS -- this is flat out wrong:
        if (SW.getSoftClipCount() > ((int) read.size() > 40 ? (int) read.size()/20 : 0))
        {
            return INT_MAX;
        }
#endif
        // slightly faster: populate the alignment vector
        SW.populateAlignment();

        if (mapperOptions.debug)
        {
            SW.debugPrint();
        }
    }
    else
    {
        SW.clearAlignment();
    }

    //
    // This is relatively speaking a rare code path, so there is
    // little benefit in tweaking it for performance - but
    // we absolutely have to have the functionality correct.
    //
    //
    // This is a little messy, so some clarification:
    // We just ran Smith Waterman backwards from before
    // the first base of the index word, and we ran
    // it after the last base, so we have three chunks
    // of information to collect - the CIGAR string for
    // the portion prior to the index word, the CIGAR string
    // for the index word, and the CIGAR string for
    // the portion following the index word.  Sometimes
    // the string before and after the index word can be
    // empty (since the index word is allowed to be at the
    // start or end of the read).
    //
    //
    //
    // this call populates the backward cigar portion from
    // the just completed H array above.
    SW.rollCigar(cigarRollerBackward);

    //
    // Here is an additional bit of messiness - the genome match
    // position that was passed in needs to be adjusted by
    // however many inserts or deletes there was (in total),
    // so that we can report the correct read match position
    // back to the user.
    //
    //
    matchPositionOffset = cigarRollerBackward.getMatchPositionOffset();

    //
    // now add in the ::wordSize matches for the index
    // word:
    //
    // NB: we don't know whether or not each base in
    // the index word matches the reference - it could
    // be that we are here due to one of the permutations.
    //
    // HOWEVER, it does not matter, the CIGAR string character
    // is 'M' regardless (M stands for match or mismatch).
    //
    cigarRollerBackward.Add(CigarRoller::match, wordIndex->wordSize);

    //
    // Now for the last portion - include the CIGAR string (ok, operations) for
    // the read:
    //
    cigarRollerBackward += cigarRollerForward;

    return 0;
}


// run a local alignment algorithm on Base Space sequences.
//
// update: localAlignmentWindowPostion
// input: localAlignmentWindowSize
// output: mismatchCount, gapOpenCount, gapCloseCount, gapExtensionCount, cigarRoller
//
genomeIndex_t ReadIndexer::localAlignment(
                                          genomeIndex_t localAlignmentWindowPostion,
                                          genomeIndex_t localAlignmentWindowSize,
                                          int &mismatchCount,
                                          int &quality,
                                          int &gapOpenCount,
                                          int &gapCloseCount,
                                          int &gapExtensionCount,
                                          CigarRoller &cigarRoller)
{
    mismatchCount = quality = gapOpenCount = gapCloseCount = gapExtensionCount = 0;

#if 0

#define MAX_SMITH_WATERMAN_REFERENCE 4096
#define MAX_SMITH_WATERMAN_READ 512

    if (localAlignmentWindowSize >= MAX_SMITH_WATERMAN_REFERENCE)
    {
        std::cerr << "FATAL: Smith Waterman Array is fixed to " << MAX_SMITH_WATERMAN_REFERENCE << "." << std::endl;
        std::cerr << "ReadIndexer::localAlignment called with localAlignmentWindowSize = " << localAlignmentWindowSize << "." << std::endl;
        std::cerr << "Either recompile with MAX_SMITH_WATERMAN_REFERENCE set to  a larger value or run Karma with a smaller value for --maxInsert" << std::endl;
        exit(1);
    }

    // XXX also check that read length is < MAX_SMITH_WATERMAN_READ

    static SmithWaterman<MAX_SMITH_WATERMAN_READ, MAX_SMITH_WATERMAN_REFERENCE, uint16_t, String , GenomeSequence , String, uint32_t, genomeIndex_t > SW;

    genomeIndex_t cigarStartingPoint = 0;
    uint32_t    softClipCount = 0;

    SW.localAlignment(
                      MAX_SMITH_WATERMAN_REFERENCE,      // XXX this is a band size parameter...
                      read,
                      read.Length(),
                      phredQuality,
                      *gs,
                      localAlignmentWindowSize,
                      localAlignmentWindowPostion,
                      cigarRoller,
                      softClipCount,
                      cigarStartingPoint,
                      quality
                      );

#if 0
    SW.debugPrint(false);
#endif

    //
    return localAlignmentWindowPostion + SW.maxCostPosition.second - SW.maxCostPosition.first;
#else
    Weight wt;
    genomeIndex_t   matchPosition = INVALID_GENOME_INDEX;
    GreedyTupleAligner< const char *, GenomeSequence &, genomeIndex_t > greedy(wt);
    greedy.Align(
                 read.c_str(),                   // query- XXXX THIS WAS A String & argument
                 read.size(),                    // query length
                 *gs,                            // reference
                 localAlignmentWindowPostion,    // index to start of reference to search
                 localAlignmentWindowSize,       // number of bases to search in reference
                 cigarRoller,                    // put the resulting cigar string here
                 matchPosition                   // put the match here
                 );
    return matchPosition + localAlignmentWindowPostion;   // plus what offsert?

#endif

}

