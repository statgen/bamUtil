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

#include "GenomeSequence.h"
#include "MappingStats.h"
#include "MapperBase.h"
#include "MapperSEColorSpace.h"
#include "ReadsProcessor.h"
#include "Error.h"
#include "MathConstant.h"
#include "Performance.h"
#include "SmithWaterman.h"
#include "Util.h"

#include "../bam/SamFlag.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

unsigned int MapperBase::colorSpaceSNP[][3]={{ 5,10,15}, { 4,11,14}, { 7, 8,13}, { 6, 9,12},
    { 1,11,14}, { 0,10,15}, { 3, 9,12}, { 2, 8,13},
    { 2, 7,13}, { 3, 6,12}, { 0, 5,15}, { 1, 4,14},
    { 3, 6, 9}, { 2, 7, 8}, { 1, 4,11}, { 0, 5,10}
};

#include "debug.h"
//
// There are several goals expressed in this matching code:
//   1) index based on the high quality end of the read,
//   2) generalize so that word 0 is always the first 'wordSize'
//      bases from the start of the original read.
//   3) generalize so that we allow a variable number of
//      'wordSize' index words, as determined by the size of
//      the read divided by 'wordSize'.
//
// Example below is two 15 character words, plus a 5 character remainder string
//
//
//   FORWARD STRAND CHECKS:
//
//                             Remainder End
//                                         |
//                                         |
//                     Remainder Start     |
//                                   |     |
//                   word1 start     |     |
//   Match Position  |               |     |
//   word0 start     |               |     |
//   |               |               |     |
//   V     word0     V     word1     V rem V
//   |_______________|_______________|_____|
//
//   Match Position  = The first base pair index in the genome sequence that
//                     corresponds to this read
//   Word0 Start     = Match Position
//   Word1 Start     = Match Position + 1 * wordSize
//   Remainder Start = Match Position + 2 * wordSize
//   Remainder End   = Match Position + read_fragment.size() - 1    (inclusive)
//
//
//
//   REVERSE STRAND CHECKS:
//
//   Match Position
//   Remainder Start
//   |
//   |     word1 start
//   |     |               word0 start
//   |     |               |
//   |     |               |       word0 end
//   |     |               |               |
//   V rem V     word1     V     word0     V
//   |_____|_______________|_______________|
//
//   Match Position  = The first base pair index in the genome sequence that
//                     corresponds to this read
//   Remainder Start = Match Position
//   Remainder End   = Match Position + Length() - 2 * wordSize - 1   (inclusive)
//   Word1 Start     = Match Position + Length() - 2 * wordSize
//   Word2 Start     = Match Position + Length() - 1 * wordSize
//
//

MatchedReadBase::~MatchedReadBase()
{
}

/*
 * update member variables using betterMatch, note, without calculating mappingquality
 * @param betterMatch a MatchedReadBase variable in which the genome match position, quality, mismatchCount, whichWord and indexer are provided
 */
void MatchedReadBase::updateMatch(MatchedReadBase &betterMatch)
{
    genomeMatchPosition = betterMatch.genomeMatchPosition;
    quality = betterMatch.quality;
    mismatchCount = betterMatch.mismatchCount;
    whichWord = betterMatch.whichWord;
    indexer = betterMatch.indexer;
}

/*
 * calcuate phred probability of mapping position is actually wrong given sumQ and posteriorProbability
 * @param sumQ              sum of phred qualities at mismatched locations
 * @posteriorProbability    given posterior probability
 */
double MatchedReadBase::getQualityScore(double sumQ, double posteriorProbability)
{
    double qs = 100.0;
    double errorProbability =
        1.0 -
        pow(10.0, (-0.1 * (sumQ))) /
        posteriorProbability;

    //
    // if errorProbability is very small, log10 of it gets quite large,
    // so top code it here:
    //
    if (errorProbability > 1e-10)
        qs = log10(errorProbability) * (-10.0);

    //
    // if errorProbability is very close to 1, qs is sometimes
    // -0,  which prints out in an ugly fashion, so bottom code it:
    //
    if (qs<1e-10) qs = 0.0;

    return qs;
}

/*
 * calcuate phred probability of mapping position is actually wrong using internal member varialbe:
 * quality and cumulativePosteriorProbabilities
 */
double MatchedReadBase::getQualityScore()
{
    if (!qualityIsValid()) return 0.0;

    return getQualityScore(quality, cumulativePosteriorProbabilities);
}


void MatchedReadBase::constructorClear()
{
    genomeMatchPosition = 0;
    quality = UNSET_QUALITY;
    mismatchCount = 0;
    cumulativePosteriorProbabilities = 0.0;
    numMatchContributors = 0;
    numBest = 0;
    whichWord = 0;
    gappedAlignment = false;
    indexer = NULL;
}

/**
 * calculate the SEQ field for base space inputs
 * @param sequce2print          stored calculation result
 * @param cigarRoller           input CigarRoller
 * @param showReferenceBases    control how to display match/mismatched bases. True: display mismatch using lower case of the original read bases; false: dispaly match using '=', others using original read bases
 * @return                      always return true
 */
bool MatchedReadBase::getBSSeq2print(std::string&    sequence2print,
                                     CigarRoller&    cigarRoller,
                                     bool            showReferenceBases)
{
    // For compression purposes, read bases that match the reference
    // are shown as '='.
    // Read bases that don't match are shown as the read version of base.
    //
    // Insertions and deletions have to be accounted for here.
    //
    // match/mismatch CIGAR operations are easy, but when
    // we get an insertion, we need to print those bases, but
    // not compare against the reference, and we need then
    // skip past the inserted bases in the read, but not skip
    // any bases in the reference.
    //
    // Similarly for deletions... Simple, no?
    //
    assert(indexer != NULL);

    if (qualityIsValid())
    {
        genomeIndex_t referencePosition = genomeMatchPosition;
        uint32_t readPosition = 0;

        //
        // Compute the SAM SEQ column.  Output is complex because it has
        // to correctly parse the CIGAR string column, iterating over each
        // CIGAR component, then over each base, printing or not printing the
        // correct output base.
        //
        // It is important to get correct, but unfortunately, it isn't clear
        // that the three or four code paths that can populate the cigarRoller
        // are always going to get it right in all circumstances.
        //
        // Therefore, we're going to add a bit of boundary checking here to terminate
        // if the cigar string appears invalid.
        //
        for (int cigarIndex = 0; cigarIndex < cigarRoller.size(); cigarIndex++)
        {
            for (uint32_t i = 0; i < cigarRoller[cigarIndex].count; i++)
            {
                char readBase = indexer->read[readPosition];
                switch (cigarRoller[cigarIndex].operation)
                {
                    case CigarRoller::match:
                    case CigarRoller::mismatch:
                    {
                        if (readBase == (*(indexer->gs))[referencePosition])
                            sequence2print.push_back(showReferenceBases ? readBase : '=');
                        else
                        {
                            // show mismatches in lower case if we're dumping all bases
                            sequence2print.push_back(showReferenceBases ? tolower(readBase) : readBase);
                        }
                        readPosition++;
                        referencePosition++;
                        break;
                    }
                    case CigarRoller::insert:
                        // for an insert, we print the base
                        sequence2print.push_back(readBase);
                        // the bases in the read are extra, so we skip over
                        // them but not the reference
                        readPosition++;
                        break;
                    case CigarRoller::del:
                        // for delete, there is nothing to print, but we skip reference bases
                        referencePosition++;
                        break;
                    case CigarRoller::softClip:
                        // for a soft clip, we always print the base
                        sequence2print.push_back(readBase);
                        readPosition++;
                        referencePosition++;
                        break;
                    case CigarRoller::hardClip:
                    case CigarRoller::pad:
                        // not implemented by CigarRoller class
                        break;
                    case CigarRoller::none:
                    default:
                        // programming error
                        assert(false);
                        break;
                }
            }
        }
    }
    else
    {
        // if the quality is invalid, trying to compare against reference is useless,
        // so at least just print the read here...
        //
        // This should be pointing to the original read (forward) by default.
        //
        sequence2print = indexer->read;
    }
    return true;
}

bool MatchedReadBase::getCSSeqAndQual2print(std::string&    sequence2print,
        std::string&    quality2print,
        std::string&    cs_read_fragment,
        std::string&    cs_data_quality,
        GenomeSequence* gs,
        GenomeSequence* csgs,
        CigarRoller&    cigarRoller,
        bool            showReferenceBases)
{
    // For compression purposes, read bases that match the reference
    // are shown as '='.
    // Read bases that don't match are shown as the read version of base.
    //
    // Insertions and deletions have to be accounted for here.
    //
    // match/mismatch CIGAR operations are easy, but when
    // we get an insertion, we need to print those bases, but
    // not compare against the reference, and we need then
    // skip past the inserted bases in the read, but not skip
    // any bases in the reference.
    //
    // Similarly for deletions... Simple, no?
    //
    assert(indexer != NULL);
    std::string read_fragment;
    std::string data_quality;


    // if showReference is true:  match will be readBase, mismatch will be lower case of readBase
    // if showReference is false: match will be  "="    , mismatch will be readBase
    if (qualityIsValid())
    {
        calibrateSequence(read_fragment, data_quality,
                          cs_read_fragment, cs_data_quality,
                          gs, csgs, genomeMatchPosition,
                          indexer->isForward);

        quality2print=data_quality;
        //        getCigarAndSequence2print(cigar, sequence2print);
        if (isForward())
        {
            for (uint32_t i = 0; i< read_fragment.size(); i++)
            {
                char readBase = read_fragment[i];
                if (readBase == (*gs)[genomeMatchPosition + i -1 ]) // correct for 0 based
                    sequence2print.push_back(showReferenceBases ? readBase : '=');
                else
                {
                    sequence2print.push_back(showReferenceBases ? tolower(readBase) : readBase);
                }
            }
        }
        else
        {
            for (uint32_t i = 0; i< read_fragment.size(); i++)
            {
                char readBase = read_fragment[i];
                if (readBase == (*gs)[genomeMatchPosition + i - 1 ])
                {
                    sequence2print.push_back(showReferenceBases ? readBase : '=');
                }
                else
                {
                    sequence2print.push_back(showReferenceBases ? tolower(readBase) : readBase);
                }
            }
            // quality2print.Invert();
        }
    }
    else
    {
        // if the quality is invalid, trying to compare against reference is useless,
        // so at least just print the read here...
        sequence2print = convertCSRead(cs_read_fragment);
        quality2print = convertCSQuality(cs_data_quality);
    }

    return true;
    if (qualityIsValid())
    {


        genomeIndex_t referencePosition = genomeMatchPosition;
        uint32_t readPosition = 0;


        //
        // Compute the SAM SEQ column.  Output is complex because it has
        // to correctly parse the CIGAR string column, iterating over each
        // CIGAR component, then over each base, printing or not printing the
        // correct output base.
        //
        // It is important to get correct, but unfortunately, it isn't clear
        // that the three or four code paths that can populate the cigarRoller
        // are always going to get it right in all circumstances.
        //
        // Therefore, we're going to add a bit of boundary checking here to terminate
        // if the cigar string appears invalid.
        //
        for (int cigarIndex = 0; cigarIndex < cigarRoller.size(); cigarIndex++)
        {
            for (uint32_t i = 0; i < cigarRoller[cigarIndex].count; i++)
            {
                char readBase = indexer->read[readPosition];
                switch (cigarRoller[cigarIndex].operation)
                {
                    case CigarRoller::match:
                    case CigarRoller::mismatch:
                    {
                        if (readBase == (*(indexer->gs))[referencePosition])
                            sequence2print.push_back(showReferenceBases ? readBase : '=');
                        else
                        {
                            // show mismatches in lower case if we're dumping all bases
                            sequence2print.push_back(showReferenceBases ? tolower(readBase) : readBase);
                        }
                        readPosition++;
                        referencePosition++;
                        break;
                    }
                    case CigarRoller::insert:
                        // for an insert, we print the base
                        sequence2print.push_back(readBase);
                        // the bases in the read are extra, so we skip over
                        // them but not the reference
                        readPosition++;
                        break;
                    case CigarRoller::del:
                        // for delete, there is nothing to print, but we skip reference bases
                        referencePosition++;
                        break;
                    case CigarRoller::softClip:
                        // for a soft clip, we always print the base
                        sequence2print.push_back(readBase);
                        readPosition++;
                        referencePosition++;
                        break;
                    case CigarRoller::hardClip:
                    case CigarRoller::pad:
                        // not implemented by CigarRoller class
                        break;
                    case CigarRoller::none:
                    default:
                        // programming error
                        assert(false);
                        break;
                }
            }
        }
    }
    else
    {
        // if the quality is invalid, trying to compare against reference is useless,
        // so at least just print the read here...
        //
        // This should be pointing to the original read (forward) by default.
        //
        sequence2print = indexer->read;
    }
    return true;
}

/**
 * calculate the NAME and POS fields in the output SAM file
 * @param chromosomeName            store calculated chromosome name
 * @param chromosomePosition        store calculated chromosome position
 * @param mateChromosomeName        store calculated mate chromosome name
 * @param mateChromosomePosition    store calculated mate chromosome position
 * @param mate                      the pointers to the other end, for paired end reads only; NULL for single end reads
 */
void MatchedReadBase::calculateNAMEandPOS(const char* &chromosomeName, int& chromosomePosition,
        const char* &mateChromosomeName, int& mateChromosomePosition,
        const MatchedReadBase* mate)
{
    // internal const variables that will be pointed to later
    static const char* defaultChromosomeName = "*";
    static const char* defaultEqualChromosomeName = "=";
    static const int defaultChromosomePosition = 0;

    // initialize results variables.
    chromosomeName = mateChromosomeName = defaultChromosomeName;
    chromosomePosition = mateChromosomePosition = defaultChromosomePosition;

    int chromosome = INVALID_CHROMOSOME_INDEX;
    int mateChromosome = INVALID_CHROMOSOME_INDEX ;
    if (qualityIsValid())
    {
        chromosome = indexer->gs->getChromosome(genomeMatchPosition);
        chromosomePosition = genomeMatchPosition - indexer->gs->getChromosomeStart(chromosome); // 0 based
        chromosomeName = indexer->gs->chromosomeName(chromosome);
    }
    if (!mate) return ;

    // deal with paired end reads
    if (mate->qualityIsValid())
    {
        mateChromosome = indexer->gs->getChromosome(mate->genomeMatchPosition);
        mateChromosomePosition = mate->genomeMatchPosition - indexer->gs->getChromosomeStart(mateChromosome); // 0 based
        mateChromosomeName = indexer->gs->chromosomeName(mateChromosome);
    }

    if (qualityIsValid() && mate->qualityIsValid())
    {
        if (mateChromosome == chromosome)
            mateChromosomeName = defaultEqualChromosomeName;
        return;
    }

    // deal with PE reads without both valid qualities
    if (!qualityIsValid() && mate && mate->qualityIsValid())
    {
        chromosomePosition = mateChromosomePosition;
        chromosomeName = mateChromosomeName;
        mateChromosomeName = defaultEqualChromosomeName;
        return;
    }
    if (qualityIsValid() && mate && !mate->qualityIsValid())
    {
        mateChromosomePosition = chromosomePosition;
        mateChromosomeName = defaultEqualChromosomeName;
        return;
    }

    // now the PE reads are both of invalid qualities
    mateChromosomeName = defaultChromosomeName;

    assert(strcmp(chromosomeName,"*") == 0);
    assert(strcmp(mateChromosomeName, "*") == 0);
    assert(chromosomePosition == 0);
    assert(mateChromosomePosition == 0);

    return ;
}

//
// Print map information in SAM format
//
// From SAM specification 0.1.0-draft:
//
// QNAME    [^ \t\n\r]+        Query pair NAME if paired; or Query NAME if unpaired 2
// FLAG    [0-9]+ [0,216-1]    bitwise FLAG (Section 2.2.2)
// RNAME    [^ \t\n\r@=]+        Reference sequence NAME 3
// POS    [0-9]+    [0,229-1]    1-based leftmost POSition/coordinate of the clipped sequence
// MAPQ    [0-9]+    [0,28-1]    MAPping Quality (phred scaled prob. that the alignment is wrong) 4
// CIGAR    ([0-9]+[MIDNSHP])+|\*        extended CIGAR string
// MRNM    [^ \t\n\r@]+        Mate Reference sequence NaMe; “=” if the same as <RNAME> 3
// MPOS    [0-9]+    [0,229-1]    1-based leftmost Mate POSition of the clipped sequence
// ISIZE    -?[0-9]+    [-229,229]    inferred Insert SIZE 5
// SEQ    [acgtnACGTN.=]+|\*        query SEQuence; “=” for a match to the reference; n/N/. for ambiguity; cases are not maintained 6,7
// QUAL    [!-~]+|\*        query QUALity; ASCII-33 gives the Phred base quality 6,7
// TAG    [A-Z][A-Z0-9]        TAG
// VTYPE    [AifZH]        Value TYPE
// VALUE    [^\t\n\r]+        match <VTYPE> (space allowed)
//
//
// 1.QNAME and FLAG are required for all alignments. If the mapping position of the query is not available, RNAME and CIGAR are set as “*”, and POS and MAPQ as 0. If the query is unpaired or pairing information is not available, MRNM equals “*”, and MPOS and ISIZE equal 0. SEQ and QUAL can both be absent, represented as a star “*”. If QUAL is not a star, it must be of the same length as SEQ.
// 2.The name of a pair/read is required to be unique in the SAM file, but one pair/read may appear multiple times in different alignment records, representing multiple or split hits. The maximum string length is 254.
// 3.If SQ is present in the header, RNAME and MRNM must appear in an SQ header record.
// 4.Providing MAPQ is recommended. If such a calculation is difficult, 255 should be applied on the assumption that the alignment is highly accurate.
// 5.If the two reads in a pair are mapped to the same reference, ISIZE equals the difference between the coordinate of the 5’-end of the mate and of the 5’-end of the current read; otherwise ISIZE equals 0. ISIZE is negative if the mate is mapped to a smaller coordinate than the current read.
// 6.Color alignments are stored as normal nucleotide alignments with additional tags describing the raw color sequences, qualities and color-specific properties (see also Note 5 in section 2.2.4).
// 7.All mapped reads are represented on the forward genomic strand. The bases are reverse complemented from the unmapped read sequence and the quality scores and cigar strings are recorded consistently with the bases. This applies to information in the mate tags (R2, Q2, S2, etc.) and any other tags that are strand sensitive. The strand field simply indicates whether this reverse complement transform was applied from the original read sequence to obtain the bases listed in the SAM file.
//
// The cigar string is passed to this method because it does not have sufficient
// information to generate it here - we need the read, read direction, match location,
// and so on, to re-run the smith waterman gapped alignment, then call the cigar roller.
//
void MatchedReadBase::print(
    std::ostream      &file,
    MatchedReadBase   *mate,
    std::string       &fragmentTag,
    bool              showReferenceBases,
    CigarRoller       &cigarRoller,
    bool              isProperAligned,
    uint16_t          samMateFlag,
    const std::string &sampleGroupID,
    const std::string &alignmentPathTag
)
{
    // check assumptions:
    assert(cigarRoller.size()>0);   // not 100% sure about this one
    assert(indexer!=NULL);
    assert(indexer->read.size() == indexer->phredQuality.size());
    if (mate)
    {
        assert(mate->indexer!=NULL);
        assert(mate->indexer->read.size() == mate->indexer->phredQuality.size());
    }

    //
    // for this read:
    //
    const char *chromosomeName = "*";
    int chromosomePosition = 0;
    // for paired reads, set the MRNM field
    // by default, we set it to "*"
    const char *mateChromosomeName = "*";
    int mateChromosomePosition = 0;

    calculateNAMEandPOS(chromosomeName, chromosomePosition,
                        mateChromosomeName, mateChromosomePosition,
                        mate);

    double mapQualityScore = 0.0;
    // if our quality is valid, map quality
    if (qualityIsValid())
    {
        mapQualityScore = getQualityScore();    // not updated after this
    }

    // Here we are printing the SEQ field of the SAM file.
    std::string sequence2print;
    getBSSeq2print(sequence2print, cigarRoller, showReferenceBases);

    // for now, we set "proper aligned flag" when
    // 1. both reads have valid quality
    // 2. both reads mapped to the same chromosomes
    isProperAligned = mate &&
                      qualityIsValid() && mate->qualityIsValid() &&
                      mateChromosomeName[0] == '=';
//
// 5.  If the two reads in a pair are mapped to the same reference, ISIZE equals the difference between the coordinate of
//     the 5ʼ-end of the mate and of the 5ʼ-end of the current read; otherwise ISIZE equals 0 (by the “5ʼ-end” we mean the
//     5ʼ-end of the original read, so for Illumina short-insert paired end reads this calculates the difference in mapping
//     coordinates of the outer edges of the original sequenced fragment). ISIZE is negative if the mate is mapped to a
//     smaller coordinate than the current read.
//

    //
    // It is pretty much up to the caller to tell us if this is an
    // aligned pair or not, since the caller knows the conditions
    // under which this can be true.
    //
    // By default, ISIZE (insertSize) is 0 - we only update it
    // here if we have a valid reason to do so.
    //
    // ISIZE is only defined (i.e. non-zero) if the reads aligned
    // to the same chromosome.
    //
    int64_t insertSize = 0;
    if (isProperAligned)
    {
        //
        // Assumption here is that both reads were a) sequenced as a
        // pair, aligned as a pair, and have valid locations
        //
        insertSize = (int64_t) genomeMatchPosition - (int64_t) mate->genomeMatchPosition;
        if (insertSize >= 0)
        {
            insertSize += indexer->read.size();
        }
        else
        {
            insertSize -= indexer->read.size();
        }

    }

    // deal with flag
    int flag = 0;
    if (qualityIsValid())
    {
        if (!isForward()) flag |= SamFlag::REVERSE;   // mark mate reverse bit
    }
    else
    {
        flag |= SamFlag::UNMAPPED; // if mate has invalid qscore, it is unmapped
    }
    if (mate)
    {
        flag |= SamFlag::PAIRED ;
        if (isProperAligned)
            flag |= SamFlag::PROPER_PAIR;     // the read is mapped in a proper pair
        flag |= samMateFlag; // mark first or second read
        if (mate->qualityIsValid())
        {
            if (!mate->isForward()) flag |= SamFlag::MATE_REVERSED;   // mark mate reverse bit
        }
        else
        {
            flag |= SamFlag::MATE_UNMAPPED; // if mate has invalid qscore, it is unmapped
        }
        // NB: we could mark flag |= 0x0200 when we get 'N' bases...
        // NB: well, not really... we can handle 'N' bases in some cases...
    }

    const char *cigarString = NULL;
    static const char* defaultCigarString = "*";

    if (qualityIsValid())
    {
        cigarString = (char *) cigarRoller.getString();
    }
    else
    {
        //
        // if this read is unmapped, we have to clear some values out
        // XXX Clean this up...
        //
        cigarString = defaultCigarString;
    }

    //
    // final sanity check:
    //  read length == qual length == cigar string M count
    //
    if (indexer->phredQuality!="*") assert(sequence2print.size() == indexer->phredQuality.size());
    if (strcmp(cigarString, "*")!=0) assert(sequence2print.size() == (unsigned int) cigarRoller.getExpectedQueryBaseCount());

    file
    << (fragmentTag.c_str() + 1)             // QNAME
    << '\t' << flag                          // FLAG
    << '\t' << chromosomeName                // RNAME
    << '\t' << chromosomePosition + 1        // POS
    << '\t' << (int)(mapQualityScore+.5)     // MAPQ (XXX needs to be %.0f)
    << '\t' << cigarString                   // CIGAR
    << '\t' << mateChromosomeName            // MRNM
    << '\t' << mateChromosomePosition + 1    // MPOS
    << '\t' << insertSize                    // ISIZE
    << '\t' << sequence2print                // SEQ
    << '\t' << indexer->phredQuality;        // QUAL

    printOptionalTags(file, isProperAligned, sampleGroupID, alignmentPathTag);

    file << std::endl;
}


//
// start = 1, end = 2
// color read:  00  (position 1, 2)
// ref genome: AACC
void MatchedReadBase::fixBaseRange(int start, int end, /// inclusive boundaries
                                   std::string& read_fragment, std::string& data_quality,
                                   const std::string& cs_read_fragment, const std::string& cs_data_quality,
                                   GenomeSequence* gs, GenomeSequence* csgs,
                                   genomeIndex_t genomeMatchPosition,
                                   bool isForwardStrand)
{
    if (genomeMatchPosition==0) return;

#define FWD_GENOME_POS(i) (genomeMatchPosition+(i)-2)
#define FWD_BSREAD_POS(i) ((i)-1)
#define BWD_GENOME_POS(i) (genomeMatchPosition+(length-1)-(i))
#define BWD_BSREAD_POS(i) ((length-1)-(i))
#define BACKWARD_CS2BS_REF_OFFSET (-1)

    // we will not fix longer (>2) color space mismatches
    // as single color space error rate is 1-2%, at least 3 consecutive error has rate 1^(-6), which is rare
    int nMismatch = end-start+1;
    if (nMismatch>2)
        return;

    // for mismatch at one position, we can infer its base call from two direction
    // one direction is from lower read index, we represent it choice0
    // the other direction is from higher read index, we represent it choice1
    char choice1, choice2;
    char qual1, qual2, qual3;
    int length = cs_read_fragment.size();

    if (isForwardStrand)
    {
        switch (nMismatch)
        {
            case 1:
                BaseAndColorToBase((*gs)[FWD_GENOME_POS(start-1)], cs_read_fragment[start], choice1);
                BaseAndColorToBase((*gs)[BWD_GENOME_POS(end)], cs_read_fragment[end], choice2);
                qual1 = data_quality[FWD_BSREAD_POS(start)];
                qual2 = data_quality[FWD_BSREAD_POS(end+1)];
                if (choice1 != choice2)
                {
                    if (qual1 > qual2)
                    {
                        read_fragment[FWD_BSREAD_POS(start)] = choice1;
                        data_quality[FWD_BSREAD_POS(start)] = qual1 - qual2 + 33;
                    }
                    else
                    {
                        read_fragment[FWD_BSREAD_POS(start)] = choice2;
                        data_quality[FWD_BSREAD_POS(start)] = qual2 - qual1 + 33;
                    }
                }
                else
                {
                    read_fragment[FWD_BSREAD_POS(start)] = choice1 ;
                    if (qual1 > qual2)
                        data_quality[FWD_BSREAD_POS(start)] =qual2;
                    else
                        data_quality[FWD_BSREAD_POS(start)] =qual1;
                }
                break;
            case 2:
                char color;
                BaseAndColorToBase((*gs)[FWD_GENOME_POS(start-1)], cs_read_fragment[start], choice1);
                BaseAndColorToBase((*gs)[FWD_GENOME_POS(end+1)], cs_read_fragment[end+1], choice2);
                BaseAndBaseToColor(choice1, choice2, color);
                if (color != cs_read_fragment[end])   // Not a SNP, so it likely there's one wrong color base.
                {
                    // we pick the two highest qualities surrounding 2 mismatch position
                    // the use these 2 colors to call bases
                    qual1 = data_quality[FWD_BSREAD_POS(start)];
                    qual2 = data_quality[FWD_BSREAD_POS(end)];
                    qual3 = data_quality[FWD_BSREAD_POS(end+1)];
                    if (qual1 >= qual3 && qual2 >= qual3)   //qual3 is the smallest
                    {
                        read_fragment[FWD_BSREAD_POS(start)] = choice1;
                        BaseAndColorToBase(choice1, cs_read_fragment[end], read_fragment[FWD_BSREAD_POS(end)]);
                        data_quality[start] = qual1 - qual3 +33;
                        data_quality[end] = qual2 - qual3 +33;
                    }
                    else if (qual1 >= qual2 && qual3 >= qual2)   //qual2 is the smallest
                    {
                        read_fragment[FWD_BSREAD_POS(start)] = choice1;
                        read_fragment[FWD_BSREAD_POS(end)] = choice2;
                        data_quality[start] = qual1 - qual2 +33;
                        data_quality[end] = qual3 - qual2 +33;
                    }
                    else   //qual1 is the smallest
                    {
                        read_fragment[FWD_BSREAD_POS(end)] = choice2;
                        BaseAndColorToBase(choice2, cs_read_fragment[end], read_fragment[FWD_BSREAD_POS(start)]);
                        data_quality[start] = qual2 - qual1 +33;
                        data_quality[end] = qual3 - qual1 +33;
                    }
                }
                else    // 2 consecutive msimatch color codes represents a SNP
                {
                    read_fragment[FWD_BSREAD_POS(start)] = choice1;
                    read_fragment[FWD_BSREAD_POS(end)] = choice2;
                    data_quality[FWD_BSREAD_POS(start)] = cs_data_quality[start];
                    data_quality[FWD_BSREAD_POS(end)] = cs_data_quality[end];
                }
                break;
            default:
                return;
        }
    }
    else   // backwardStrand
    {
        switch (nMismatch)
        {
            case 1:
                BaseAndColorToBase((*gs)[BWD_GENOME_POS(start-1)], cs_read_fragment[start], choice1);
                BaseAndColorToBase((*gs)[BWD_GENOME_POS(end)], cs_read_fragment[end], choice2);
                qual1 = data_quality[BWD_BSREAD_POS(start)];
                qual2 = data_quality[BWD_BSREAD_POS(end+1)];
                if (choice1 != choice2)
                {
                    if (qual1 > qual2)
                    {
                        read_fragment[BWD_BSREAD_POS(start)] = choice1;
                        data_quality[BWD_BSREAD_POS(start)] = qual1 - qual2 + 33;
                    }
                    else
                    {
                        read_fragment[BWD_BSREAD_POS(start)] = choice2;
                        data_quality[BWD_BSREAD_POS(start)] = qual2 - qual1 + 33;
                    }
                }
                else
                {
                    read_fragment[BWD_BSREAD_POS(start)] = choice1 ;
                    if (qual1 > qual2)
                        data_quality[BWD_BSREAD_POS(start)] =qual2;
                    else
                        data_quality[BWD_BSREAD_POS(start)] =qual1;
                }
                break;
            case 2:
                char color;
                BaseAndColorToBase((*gs)[BWD_GENOME_POS(start-1) + BACKWARD_CS2BS_REF_OFFSET],
                                   cs_read_fragment[start], choice1);
                BaseAndColorToBase((*gs)[BWD_GENOME_POS(end+1) +BACKWARD_CS2BS_REF_OFFSET],
                                   cs_read_fragment[end+1], choice2);
                BaseAndBaseToColor(choice1, choice2, color);
                if (color != cs_read_fragment[end])
                {
                    // we pick the two highest qualities surrounding 2 mismatch position
                    // the use these 2 colors to call bases
                    qual1 = data_quality[BWD_BSREAD_POS(start)];
                    qual2 = data_quality[BWD_BSREAD_POS(end)];
                    qual3 = data_quality[BWD_BSREAD_POS(end+1)];
                    if (qual1 >= qual3 && qual2 >= qual3)   //qual3 is the smallest
                    {
                        read_fragment[BWD_BSREAD_POS(start)] = choice1;
                        BaseAndColorToBase(choice1, cs_read_fragment[end], read_fragment[BWD_BSREAD_POS(end)]);
                        data_quality[start] = qual1 - qual3 +33;
                        data_quality[end] = qual2 - qual3 +33;
                    }
                    else if (qual1 >= qual2 && qual3 >= qual2)   //qual2 is the smallest
                    {
                        read_fragment[BWD_BSREAD_POS(start)] = choice1;
                        read_fragment[BWD_BSREAD_POS(end)] = choice2;
                        data_quality[start] = qual1 - qual2 +33;
                        data_quality[end] = qual3 - qual2 +33;
                    }
                    else   //qual1 is the smallest
                    {
                        read_fragment[BWD_BSREAD_POS(end)] = choice2;
                        BaseAndColorToBase(choice2, cs_read_fragment[end], read_fragment[BWD_BSREAD_POS(start)]);
                        data_quality[start] = qual2 - qual1 +33;
                        data_quality[end] = qual3 - qual1 +33;
                    }
                }
                else    // 2 consecutive msimatch color codes represents a SNP
                {
                    read_fragment[BWD_BSREAD_POS(start)] = choice1;
                    read_fragment[BWD_BSREAD_POS(end)] = choice2;
                    data_quality[BWD_BSREAD_POS(start)] = cs_data_quality[start];
                    data_quality[BWD_BSREAD_POS(end)] = cs_data_quality[end];
                }
                break;
            default:
                return;
        }
    }

}


/**
 * When translate color space reads to base space, it is possible we have long strands of mismatch
 * therefore, we will try to solve this problem here:
 * mark the start and end position for the consecutive mismatches,
 * then call fixBaseRange() function.
 * @param read_fragment         store translated read in base space
 * @param data_quality          store translated quality in base space
 * @param cs_read_fragment      original read in color space
 * @param cs_data_quality       original quality in color space
 * @param gs                    a pointer to the base space reference genome
 * @param csgs                  a pointer to the color space reference genome
 * @param genomeMatchPosition   the match position
 * @param isForwardStrang       strand direction
 */
void MatchedReadBase::calibrateSequence(std::string& read_fragment, std::string& data_quality,
                                        const std::string& cs_read_fragment, const std::string& cs_data_quality,
                                        GenomeSequence* gs, GenomeSequence* csgs,
                                        genomeIndex_t genomeMatchPosition,
                                        bool isForwardStrand)
{
    int length = cs_read_fragment.size();
    int cumMismatch = 0;
    int cumMismatchStart = 0;
    int cumMismatchEnd = 0;

//    read_fragment = String('N',length-1);
//    read_fragment = String('N',length-1).c_str();
    read_fragment.resize(length-1);
    data_quality.resize(length-1);
    fill(read_fragment.begin(), read_fragment.begin()+length-1, 'N');
    if (isForwardStrand)
        data_quality = cs_data_quality.substr(1);
//        data_quality = cs_data_quality.Right(length-1);
    else
    {
        std::string tmp = cs_data_quality.substr(1);
        reverse(tmp.begin(), tmp.end());
        data_quality = tmp;
//        data_quality = cs_data_quality.Right(length-1).Invert();
    }

    if (genomeMatchPosition == INVALID_GENOME_INDEX)
        return;

    // Deal with the base next to primer.
    //
    // XXX the problem is that the primer base means nothing, according to
    // Tom, and I think at least one online source I saw.
    //

    // We always set the base next to primer to be a colorspace-basespace translation,
    // and assume it is the correct base call.
    if (isForwardStrand)
        BaseAndColorToBase(cs_read_fragment[0], cs_read_fragment[1], read_fragment[0]);
    else
    {
        BaseAndColorToBase(cs_read_fragment[0], cs_read_fragment[1], read_fragment[length-1-1]);
        read_fragment[length-1-1] = GenomeSequence::base2complement[
                                        (int)(read_fragment[length-1-1])];
    }
    // Notice for clarification of indices:
    ////////////////////////////////////////////////////
    // for a forward strand match,
    // for index i: 0 1 2 3 4 5 ...
    // color read:  A 0 1 2 3 0 ...
    // color ref:       1 2 3 0 ...
    //                  ^       genomeMatchPosition + i-2
    // bs read:       A C G T T ...
    //                  ^       read_fragment[i-1]

    // for a backward strand match,
    // for index i: (len-1).. 5 4 3 2 1 0
    // color read:  0  ...... 0 1 2 3 0 A
    // color ref:   0  ...... 0 1 2 3
    //              ^       genomeMatchPosition + (len-1) - i
    // bs read:     .....     C C G A A
    //              ^       read_fragment[(len-1)-i]

    char csRefBase;
    char bsRefBase;
    for (int i = 2; i < length; ++i)
    {
        if (isForwardStrand)
            csRefBase = (*csgs)[FWD_GENOME_POS(i)];
        else
            csRefBase = (*csgs)[BWD_GENOME_POS(i)];
        if (cs_read_fragment[i] != csRefBase)
        {
            if (cumMismatch == 0)
            {
                cumMismatch ++;
                cumMismatchStart = i;
                cumMismatchEnd = i;
            }
            else
            {
                cumMismatch ++;
                cumMismatchEnd ++;
            }
        }
        else
        {
            if (isForwardStrand)
            {
                bsRefBase = (*gs)[FWD_GENOME_POS(i)];
                read_fragment[i-1] = bsRefBase;
            }
            else
            {
                bsRefBase = (*gs)[BWD_GENOME_POS(i) + BACKWARD_CS2BS_REF_OFFSET];
                read_fragment[(length-1)-i] = bsRefBase;
            }

            if (cumMismatch != 0)
            {
                cumMismatch = 0;
                fixBaseRange(cumMismatchStart, cumMismatchEnd,
                             read_fragment, data_quality,
                             cs_read_fragment, cs_data_quality,
                             gs, csgs, genomeMatchPosition,
                             isForwardStrand);
            }
        }
    }

    if (cumMismatch != 0)
    {
        // we manually fix the last base call and
        // set its quality to be empirically half of its value
        if (isForwardStrand)
        {
            bsRefBase = (*gs)[FWD_GENOME_POS(length-1)];
            read_fragment[FWD_BSREAD_POS(length-1)] = bsRefBase;
            data_quality[(length-1)-1] = (data_quality[(length-1)-1]-33)/2 + 33;
        }
        else
        {
            bsRefBase = (*gs)[BWD_GENOME_POS(length-1)];
            read_fragment[0] = bsRefBase;
            data_quality[0] = (data_quality[0]-33)/2 +33;
        }
        fixBaseRange(cumMismatchStart, length-2, //length-2: the second last base of cs_read_fragment
                     read_fragment, data_quality,
                     cs_read_fragment, cs_data_quality,
                     gs, csgs, genomeMatchPosition,
                     isForwardStrand);
    }
}


void MatchedReadBase::printColorSpace(
    std::ostream &file,
    GenomeSequence* gs,
    GenomeSequence* csgs,
    MatchedReadBase *mate,
    std::string &cs_read_fragment,
    std::string &cs_data_quality,
    std::string &fragmentTag,
    bool showReferenceBases,
    CigarRoller  &cigarRoller,
    bool        isProperAligned,
    uint16_t        samMateFlag,
    const std::string     &sampleGroupID,
    const std::string     &alignmentPathTag
)
{
    assert(!gs->isColorSpace());
    assert(csgs->isColorSpace());  // color space reads

#if 0
    // here we insert debug code for testing!!!
    genomeMatchPosition = 1225662037;
    cs_read_fragment = "A30330323211222230003322012223031233";
    //_read_fragment = "A30330323233222230003322012223031233";
    cs_data_quality =   "!1111111111111111111111111111111111";
    // end test code
#endif
    // check assumptions:
    // std::cout << "cigarRoller.getString()=" << cigarRoller.getString() << std::endl;
    assert(cigarRoller.size()>0);   // not 100% sure about this one
    assert(indexer!=NULL);
    assert(indexer->read.size() == indexer->phredQuality.size());
    if (mate)
    {
        assert(mate->indexer!=NULL);
        assert(mate->indexer->read.size() == mate->indexer->phredQuality.size());
    }

    //
    // for this read:
    //
    const char *chromosomeName = "*";
    int chromosomePosition = 0;
    // for paired reads, set the MRNM field
    // by default, we set it to "*"
    const char *mateChromosomeName = "*";
    int mateChromosomePosition = 0;

    calculateNAMEandPOS(chromosomeName, chromosomePosition,
                        mateChromosomeName, mateChromosomePosition,
                        mate);

    double mapQualityScore = 0.0;
    // if our quality is valid, map quality
    if (qualityIsValid())
    {
        mapQualityScore = getQualityScore();    // not updated after this
    }

    // Here we are printing the SEQ field of the SAM file.
    std::string sequence2print;
    std::string quality2print;
    D_STRING(cs_read_fragment);
    getCSSeqAndQual2print(sequence2print, quality2print,
                          cs_read_fragment, cs_data_quality,
                          gs, csgs,
                          cigarRoller, showReferenceBases);

    D_STRING(sequence2print);
    // for now, we set "proper aligned flag" when
    // 1. both reads have valid quality
    // 2. both reads mapped to the same chromosomes
    isProperAligned = mate &&
                      qualityIsValid() && mate->qualityIsValid() &&
                      mateChromosomeName[0] == '=';

    int64_t insertSize = 0;
    if (isProperAligned)
    {
        //
        // Assumption here is that both reads were a) sequenced as a
        // pair, aligned as a pair, and have valid locations
        //
        insertSize = (int64_t) genomeMatchPosition - (int64_t) mate->genomeMatchPosition;
        if (insertSize >= 0)
        {
            insertSize += indexer->read.size();
        }
        else
        {
            insertSize -= indexer->read.size();
        }

    }

    // deal with flag
    int flag = 0;
    if (qualityIsValid())
    {
        if (!isForward()) flag |= SamFlag::REVERSE;   // mark mate reverse bit
    }
    else
    {
        flag |= SamFlag::UNMAPPED; // if mate has invalid qscore, it is unmapped
    }
    if (mate)
    {
        flag |= SamFlag::PAIRED ;
        if (isProperAligned)
            flag |= SamFlag::PROPER_PAIR;     // the read is mapped in a proper pair
        flag |= samMateFlag; // mark first or second read
        if (mate->qualityIsValid())
        {
            if (!mate->isForward()) flag |= SamFlag::MATE_REVERSED;   // mark mate reverse bit
        }
        else
        {
            flag |= SamFlag::MATE_UNMAPPED; // if mate has invalid qscore, it is unmapped
        }
    }

    const char *cigarString = NULL;
    static const char* defaultCigarString = "*";

    if (qualityIsValid())
    {
        // TODO
        // now manually add '1M' to the left of the current cigar roller, is there better way?
        CigarRoller temp;
        temp.Add(CigarRoller::match,1);
        temp.Add(cigarRoller);
        cigarString = (char *) temp.getString();
    }
    else
    {
        cigarString = defaultCigarString;
    }

    //
    // final sanity check:
    //  read length == qual length == cigar string M count
    if (indexer->phredQuality!="*") assert(sequence2print.size() == quality2print.size());
    if (strcmp(cigarString, "*")!=0) assert(sequence2print.size() == 1 + (unsigned int) cigarRoller.getExpectedQueryBaseCount());

    file
    << (fragmentTag.c_str() + 1)             // QNAME
    << '\t' << flag                          // FLAG
    << '\t' << chromosomeName                // RNAME
    << '\t' << chromosomePosition /*+ 1 */   // POS ( here -1 is for adjust the difference between CS ref and BS ref
    << '\t' << (int)(mapQualityScore+.5)     // MAPQ (XXX needs to be %.0f)
    << '\t' << cigarString                   // CIGAR
    << '\t' << mateChromosomeName            // MRNM
    << '\t' << mateChromosomePosition + 1    // MPOS
    << '\t' << insertSize                    // ISIZE
    << '\t' << sequence2print                // SEQ
    << '\t' << quality2print;                // QUAL

    printOptionalTags(file, isProperAligned, sampleGroupID, alignmentPathTag);

    file << "\tCS:Z:"<<cs_read_fragment;
    file << "\tCQ:Z:"<<cs_data_quality;

    file << std::endl;
}

/**
 * If trimming happened, restore to original read and quality and adjust genomeMatchingPosition
 * @param cigarRoller   store cigarRoller results
 * @param indexer       change read and quality in the indexer determined by bestMatch
 * @param genomeMatchPosition   store adjusted genomeMatchPosition
 */
void MapperBase::restoreTrimming(CigarRoller&   cigarRoller,
                                 ReadIndexer&   indexer,
                                 genomeIndex_t& genomeMatchPosition)
{
    // if trimming never occured, just return
    if (!(indexer.leftTrimmedBases.size() || indexer.rightTrimmedBases.size()))
        return;

    // rewrite the cigar string to include the two softclips
    CigarRoller newCigar;

    // recover the untrimmed read and quality string:
    std::string originalRead;
    std::string originalQualities;

    if (indexer.isForward)
    {
        originalRead = indexer.leftTrimmedBases + indexer.read + indexer.rightTrimmedBases;
        originalQualities = indexer.leftTrimmedQualities + indexer.phredQuality + indexer.rightTrimmedQualities;
        newCigar.Add(CigarRoller::softClip, indexer.leftTrimmedBases.size());
        newCigar.Add(cigarRoller);
        newCigar.Add(CigarRoller::softClip, indexer.rightTrimmedBases.size());
        // conditionally adjust the match postion for the soft clip on the left
        if (genomeMatchPosition != INVALID_GENOME_INDEX)
            genomeMatchPosition -= indexer.leftTrimmedBases.size();

    }
    else
    {
        originalRead = indexer.rightTrimmedBases + indexer.read + indexer.leftTrimmedBases;
        originalQualities = indexer.rightTrimmedQualities + indexer.phredQuality + indexer.leftTrimmedQualities;
        newCigar.Add(CigarRoller::softClip, indexer.rightTrimmedBases.size());
        newCigar.Add(cigarRoller);
        newCigar.Add(CigarRoller::softClip, indexer.leftTrimmedBases.size());
        // conditionally adjust the match postion for the soft clip on the left
        if (genomeMatchPosition != INVALID_GENOME_INDEX)
            genomeMatchPosition -= indexer.rightTrimmedBases.size();
    }

    cigarRoller = newCigar;
    indexer.read = originalRead;
    indexer.phredQuality = originalQualities;
}



//
// get the CIGAR string for this match position.
//
// Returns a match position offset that indicates how to adjust the
// given match position.  This allows for variations due to indels
// that may occur ahead of where the index word position is.
//
// When the match being printed (usually MapperBase::bestMatch) was
// a gapped match (i.e. Smith Waterman was run to find it), then we
// have a messy problem.
//
// The difficulty is that to generate the CIGAR string, we have to have
// the array H, which gets discarded after computing getSumQSW().  The
// reason is we keep possibly thousands of copies of the MatchedReadBase
// in paired end code, and the array is large and is statically allocated,
// so we can't expect it to be kept around very long.
//
// In order to obtain the CIGAR string, we have to rerun the
// Smith Waterman algorithm so we can do the walkback from the lower
// right portion of the array H.
//
//
void  MapperBase::populateCigarRollerAndGenomeMatchPosition()
{
    MatchedReadBase &bestMatch = getBestMatch();    // getBestMatch is virtual, remember

    if (bestMatch.qualityIsValid() && bestMatch.gappedAlignment)
    {
        int matchPositionOffset = 0;

        //
        // if gapped, and if a valid quality, go ahead and get the
        // CIGAR roller data structure.  If quality is invalid, the code
        // below is not safe to call.
        //
        // This code fixes the issue that the index alignment will
        // fix a read at a certain location, but in the event of a
        // structural variation that occurs before the index word used
        // to assign gemomeMatchPosition for the read, we need to
        // adjust the real read location by the amount of that
        // structural variation (zero or more, actually).
        //
        // In the calls below, bestMatch.whichWord tells us which
        // index word was used to assign genomeMatchPosition, and
        // therefore where we need to adjust to.
        //
        // For a visualization of this, imagine this read:
        //
        // ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG
        //        ^       ^
        //        |       this is an example location where the index word starts (base #15)
        //        +- we'll imagine there is a single base deletion here
        // On fast index alignment, the index code places the start
        // of the read one single base earlier than it should be.
        //
        // Only by examining the cigar string can we know what the
        // adjustment is to be made.
        //
        // It is true that the Smith Waterman algorithm could incorporate
        // an offset argument that it would update according to its detection
        // of indels, but this adds even more complexity to an already complex
        // algorithm.
        //
        if (bestMatch.isForward())
        {
            //
            // To get the cigar roller, the following call is
            // essentially re-running Smith Waterman again, exactly
            // as it was done during the fast index lookup phase.
            //
            forward.getCigarRoller(bestMatch.genomeMatchPosition,
                                   bestMatch.whichWord,
                                   matchPositionOffset,
                                   cigarRoller);
        }
        else
        {
            backward.getCigarRoller(bestMatch.genomeMatchPosition,
                                    bestMatch.whichWord,
                                    matchPositionOffset,
                                    cigarRoller);
        }
        bestMatch.genomeMatchPosition += matchPositionOffset;
    }
    else if (!localGappedAlignment)
    {
        //
        // ungapped, or no match case, so don't update match position:
        //
        // Set a simple CIGAR string, basically "%dM" where %d
        // is the read length.
        //
        cigarRoller.clear();
        cigarRoller.Add(CigarRoller::match, forward.read.size());
    }
    else
    {
        // local gapped alignment yields a correct CIGAR string
        // and match position, so nothing to do here.
    }

    //
    // now unwind quality clipping changes that we made in ReadInder::setReadAndQuality
    //
    // re-assemble the indexer read (sequence) and quality strings and
    // also adjust genomeMatchPostion and the cigarRoller...
    //

    if (!bestMatch.qualityIsValid())
    {
        // XXX this is a hack to fix an earlier error
        // to fix it, put an assert() here and find out
        // what is causing the bug.
        bestMatch.indexer = &forward;
    }

    ReadIndexer &indexer = bestMatch.isForward() ? forward : backward;
    restoreTrimming(cigarRoller, indexer, bestMatch.genomeMatchPosition);
}

void MatchedReadBase::printOptionalTags(std::ostream &file, bool isPairAligned, const std::string &sampleGroupID, const std::string &alignmentPathTag)
{
    if (alignmentPathTag != "") file << "\tXA:Z:" << alignmentPathTag;

    if (sampleGroupID != "") file << "\tRG:Z:" << sampleGroupID;

    if (numMatchContributors>1) file << "\tHA:i:" << numMatchContributors;  // HA -> Hits, All

    if (qualityIsValid() && quality!=0) file << "\tUQ:i:" << quality;       // phred quality of this read, assuming it is mapped correctly

    if (numBest>1) file << "\tNB:i:" << numBest;

    //
    // XXX this should be optional:
    //
    if (!qualityIsValid())
    {
        switch (quality)
        {
            case UNSET_QUALITY:
                file << "\tER:Z:no_match";
                break;
            case INVALID_DATA:
                file << "\tER:Z:invalid_bases";
                break;
            case EARLYSTOP_DUPLICATES:
                file << "\tER:Z:duplicates";
                break;
            case EARLYSTOP_QUALITY:
                file << "\tER:Z:quality";
                break;
            case REPEAT_QUALITY:        // should never occur now
                file << "\tER:Z:repeats";
                break;
        }
    }
//    fprintf(file, "\tGP:i:%u", genomeMatchPosition);        // for debugging
}

void MatchedReadBase::debugPrint()
{
    cout << (isForward() ? 'F' : 'R')
         << "\t"
         << whichWord
         << "\t"
         << quality
         << "\t"
         << mismatchCount
         << "\t"
         << genomeMatchPosition
         << std::endl;
}

MapperBase::MapperBase() :  forward(mapperOptions), backward(mapperOptions)
{
    for (int i = 0; i < MAX_Q2P; i++)
        sumQualityToProb[i] = pow(10.0, (-0.1 * i));
    rand = &globalRandom;
    gs = NULL;
    wordIndex = NULL;
    line = 0;
    numberReads = 0;

    // lame:
    forward.isForward = true;
    backward.isForward = false;

    localGappedAlignment = false;
    isProperAligned = false;

    samMateFlag = 0;        // for paired end, mark which is mate 1 versus mate 2
}

MapperBase::~MapperBase()
{
    // we don't own gs or wordIndex, so we don't close them here.
}

int MapperBase::Word2Integer(std::string & word, unsigned int index, int &countNbases)
{
    assert(word.size() + index >= wordIndex->wordSize && index >= 0);

    int crtbaseinteger;
    int wordinteger = 0;

    //
    // XXX this is classic indication of a need to subclass, but I don't
    // know how practical it is for karma.
    //
    if (gs->isColorSpace())
    {
        //
        // NB: the + 1 is because the first character of the
        // color space read is a base pair key which we ignore here.
        // This is not a clean solution ... find something better.
        //
        for (unsigned int i = index + 1; i < index + wordIndex->wordSize + 1; i++)
        {
            wordinteger <<= 2;
            if (!isdigit(word[i])) return INVALID_WORDINDEX;
            crtbaseinteger = word[i] & 0xf; // works for ASCII and EBCDIC! whoo!
            wordinteger |= crtbaseinteger;
        }
    }
    else
    {
        for (unsigned int i = index; i < index + wordIndex->wordSize; i++)
        {
            wordinteger <<= 2;
            crtbaseinteger = GenomeSequence::base2int[(int) word[i]];

            switch (crtbaseinteger)
            {
                case GenomeSequence::baseXIndex:
                    return INVALID_WORDINDEX;   // policy issue - should a read with a single invalid base be tossed?
                case GenomeSequence::baseNIndex:
                    crtbaseinteger = 0;         // knowing this is 0, we'll visit all edits of this base later
                default:
                    wordinteger |= crtbaseinteger;
            }
        }
    }

    return wordinteger;
}

std::string MapperBase::Integer2Word(wordInteger_t n, unsigned int wordsize)
{
    assert(wordsize == wordIndex->wordSize);

    std::string word("NNNNN00000NNNNN");

    if (n > THIRTY_BIT_MASK)
        error("Invalid integer for a word of size %u: %d\n", wordsize, n);

    for (unsigned int i = 0; i < wordsize; i++)
    {
        word[wordsize - 1 - i] = MapIntegerToBase(n & 3);
        n >>= 2;
    }

    return word;
}

bool MapperBase::setReadAndQuality(const char *r, int len, const char *q)
{

    originalRead = r;
    originalQuality = q;

    // for color space reads,
    // we will convert char* r to base pair space
    // notice, color space and base space are treated differently:
    // base space:  forwardRead="ACGTA" backwardRead="ACGTA"
    // color space: forwardRead="ACGTA" backwardRead="TGCAT"
    // since if color space read is 01230, then reverse complement is 03210,
    // after map 0->A, ..., 3->T, actually, we are mapping forward ACGTA, the reverse complement ATGCA,
    // so we set the forward read to ACGTA, and set backward read to ATGCA.


    if (gs->isColorSpace() && isalpha(r[0]))
    {
        // Clip the primer base - it is not useful in our mapping code.
        // also clip the first color as it represents transition information
        // between the primer base and the first real base.
        forward.setReadAndQuality(r + 2, len - 2, q + 2);
        backward.setReadAndQuality(r + 2, len - 2, q + 2);
    }
    else
    {
        // base space - no translations are necessary - setReadAndQuality takes
        // care of converting the read to the backwards strand
        forward.setReadAndQuality(r,len,q);
        backward.setReadAndQuality(r,len,q);
    }

    //
    // Here we check the return code - if true, it means
    // we have too few valid index words for this read.
    //
    if (forward.setIndexStrategy() || backward.setIndexStrategy())
    {
        // when failing, make sure we clear out old data.
        getBestMatch().constructorClear();
        getBestMatch().indexer = &forward;
        return true;
    }
    return false;
}


//
// Given the short read input file in FASTQ format, populate the read fragment,
// data quality, fragment title (header/tag/whatever), and update the line #
// and # of reads.  Consider refactoring - a little messy looking.
//
// Returns:
//   1->read too short (deprecated)
//   2->read and quality lengths differ
//   3->too few index words
//
// Aborts on mal-formed FASTQ file.  This is a bit harsh, but does at least
// prevent ugly termination later when trying to read, e.g. a binary file.
//
int MapperBase::getReadAndQuality(IFILE f)
{
    std::string ignore;
    std::string dataQuality;
    std::string readFragment;

    while (!ifeof(f))
    {
        line ++;

        // read line for readFragment name;
        f >> fragmentTag;

        // looking for first non-empty line
        if (fragmentTag.size() > 0) break;
    }

    if (ifeof(f)) return EOF;

    numberReads ++;

    if (fragmentTag[0] != '@')
        error("Tag line must start with [@]. check line %lu\n", line);

    // read actual readFragment
    line ++;
    f >> readFragment;
    if (readFragment.size() == 0)
        error("unexpected empty short read DNA strand line %lu\n", line);

    line++;
    f >> ignore; // ignore - it's a repeat of first line

    line++;
    f >> dataQuality; // quality data
    if (dataQuality.size() == 0)
        error("unexpected empty DNA quality line %lu\n", line);

#if 0
    // debug_by_tag
    // for debug a certain read
    String tag = fragmentTag.c_str();
    String toFind = "1028:13165";
    if (tag.Find(toFind) >= 0)
    {
        printf("1");
    }
#endif

    return processReadAndQuality(fragmentTag, readFragment, dataQuality);
}

//
// Given string representations of the read, set up all internal data structures
// neccesary to perform mapping.
//
int MapperBase::processReadAndQuality(std::string& fragmentTag, std::string& readFragment, std::string& dataQuality)
{
    this->fragmentTag=fragmentTag;
    //
    // left and right truncate if requested
    //
    // It is possible that after trimming, the read
    // will be too short to map, but it is the job of
    // ::setReadAndQuality to determine that.
    //
    if (mapperOptions.trimLeft)
    {
        readFragment.erase(0, mapperOptions.trimLeft);
        dataQuality.erase(0, mapperOptions.trimLeft);
    }
    if (mapperOptions.trimRight)
    {
        throw std::logic_error("trim right option not tested");
        readFragment.erase(readFragment.size() - mapperOptions.trimRight);
        dataQuality.erase(dataQuality.size() - mapperOptions.trimRight);
    }

    // keep checks simple for now:
    if (readFragment.size() != dataQuality.size())
    {
        return 2;
    }

    //
    // on output of colorspace reads, we need the original read
    // for output to the SAM file.
    //
    if (gs->isColorSpace())
    {
        // save the read with the primer base
        ((MapperSEColorSpace *)this)->originalCSRead=readFragment;
        ((MapperSEColorSpace *)this)->originalCSQual=dataQuality;
    }

    // let setReadAndQuality conditionally truncate the first primer base
    // for color space reads
    //
    // setReadAndQuality has to be called once per read before
    // anything else is done with the class
    //

    if (setReadAndQuality(readFragment.c_str(), readFragment.size(), dataQuality.c_str()))
        return 1;

    return 0;
}

void MapperBase::debugPrint(MatchedReadBase &matchedRead)
{
    std::string dataQuality;

    for (vector<uint8_t>::iterator it = forward.binaryQuality.begin() ; it < forward.binaryQuality.end(); it++)
        dataQuality.push_back(*it);

    if (!matchedRead.isForward())
        reverse(dataQuality.begin(), dataQuality.end());

    gs->debugPrintReadValidation(
        matchedRead.indexer->read,
        dataQuality,
        matchedRead.isForward() ? 'F' : 'R',
        matchedRead.genomeMatchPosition,
        matchedRead.quality,
        matchedRead.mismatchCount,
        false);
}

//
// XXX we need two distinct modes of operation.
//
// 1) recurse to evaluate a single position
// 2) recurse to evaluate a group of positions
//
// We could:
// 1) add a second argument, evalFunctionType2, and if non-zero,
// we call it instead.
// 2) we pass a void * and a boolean to determine which one it is
// 3) turn these all into macros, which makes debug a nightmare,
//    but does solve the problem.
//
#define CANDIDATE_LIMIT 500
#if defined(CANDIDATE_LIMIT)
static int totalCandidateCount;
#endif

// a wrapper function
// for calling "evalBaseSpaceRead()" for forward & backward ReadIndex
void MapperBase::evalBaseSpaceReads(
    evalSinglePositionFunctionType onePositionMethod,
    evalAllPositionsFunctionType allPositionsMethod
)
{
    assert((onePositionMethod==NULL) ^(allPositionsMethod==NULL));
#if defined(CANDIDATE_LIMIT)
    totalCandidateCount = 0;
#endif
    if (evalBaseSpaceRead(onePositionMethod, allPositionsMethod, forward)) return;
    if (evalBaseSpaceRead(onePositionMethod, allPositionsMethod, backward)) return;
}

// a wrapper function
// for calling "evalAllCandidatesForWord()" function on
// 1. every word stored in ReadIndexer class
// 2. every word stored (replace "N" with a letter) in ReadIndexer class
// 3. (optional) mutated every word stored in ReadIndexer class
bool MapperBase::evalBaseSpaceRead(
    evalSinglePositionFunctionType onePositionMethod,
    evalAllPositionsFunctionType allPositionsMethod,
    ReadIndexer &indexer)
{
    for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
    {
        if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, 0)) return true;
        if (indexer.wordNLocations[whichWord] > 0)
        {
            wordInteger_t mask1 = 1 << indexer.wordNLocations[whichWord];
            wordInteger_t mask2 = 2 << indexer.wordNLocations[whichWord];
            wordInteger_t mask3 = 3 << indexer.wordNLocations[whichWord];

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
        }
    }

    if (indexer.checkEdits)
    {

        wordInteger_t mask1 = 1;
        wordInteger_t mask2 = 2;
        wordInteger_t mask3 = 3;
        for (unsigned int i=0; i<wordIndex->wordSize; mask1<<=2, mask2<<=2, mask3<<=2, i++)
        {
            for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
            {

                // reads may now have invalid index words, but we
                // won't use them because of course they are invalid.
                //

                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;

                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;

                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
            }
        }
    }
    return false;   // keep going
}

// a wrapper function
// for calling "evalBaseSpaceRead()" for forward & backward ReadIndex
void MapperBase::evalColorSpaceReads(
    evalSinglePositionFunctionType onePositionMethod,
    evalAllPositionsFunctionType allPositionsMethod
)
{
    assert((onePositionMethod!=NULL) || (allPositionsMethod!=NULL));
    if (evalColorSpaceRead(onePositionMethod, allPositionsMethod, forward)) return;
    if (evalColorSpaceRead(onePositionMethod, allPositionsMethod, backward)) return;
}

// a wrapper function
// for calling "evalAllCandidatesForWord()" function on
// 1. every word stored in ReadIndexer class
// 2. every word stored (replace "N" with a letter) in ReadIndexer class
// 3. (optional) mutated every word stored in ReadIndexer class
bool MapperBase::evalColorSpaceRead(
    evalSinglePositionFunctionType onePositionMethod,
    evalAllPositionsFunctionType allPositionsMethod,
    ReadIndexer &indexer)
{
    for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
    {
        if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, 0)) return true;
        if (indexer.wordNLocations[whichWord] > 0)
        {
            wordInteger_t mask1 = 1 << indexer.wordNLocations[whichWord];
            wordInteger_t mask2 = 2 << indexer.wordNLocations[whichWord];
            wordInteger_t mask3 = 3 << indexer.wordNLocations[whichWord];

            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;
            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;
            if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
        }
    }

    if (indexer.checkEdits)
    {
        wordInteger_t mask1 = 1;
        wordInteger_t mask2 = 2;
        wordInteger_t mask3 = 3;
        for (unsigned int i=0; i<wordIndex->wordSize; mask1<<=2, mask2<<=2, mask3<<=2, i++)
        {
            for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
            {

                // reads may now have invalid index words, but we
                // won't use them because of course they are invalid.
                //
                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask1)) return true;
                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask2)) return true;
                if (evalAllCandidatesForWord(onePositionMethod, allPositionsMethod, indexer, whichWord, mask3)) return true;
            }
        }

        // check single site SNPs
        unsigned int shift = 0;
        for (unsigned int i=0; i< wordIndex->wordSize; i++, shift+=2)
        {
            for (unsigned int whichWord=0; whichWord<indexer.wordInts.size(); whichWord++)
            {
                if (evalAllCandidatesForColorSpaceWord(onePositionMethod, allPositionsMethod, indexer, whichWord, shift, 0)) return true;
                if (evalAllCandidatesForColorSpaceWord(onePositionMethod, allPositionsMethod, indexer, whichWord, shift, 1)) return true;
                if (evalAllCandidatesForColorSpaceWord(onePositionMethod, allPositionsMethod, indexer, whichWord, shift, 2)) return true;
            }
        }
    }
    return false;   // keep going
}

// a wrapper function
// to call either "allPositionMethod" (a function pointer passed in as parameter)
// or calling "evalAllCandidatePositions" function using parameter "onePositionMethod"
bool MapperBase::evalAllCandidatesForWord(
    evalSinglePositionFunctionType onePositionMethod,
    evalAllPositionsFunctionType allPositionsMethod,
    ReadIndexer &indexer,
    unsigned int whichWord,
    wordInteger_t xorMask)
{
    int count = 0;
    genomeIndex_t *candidates = NULL;

    if (wordIndex->wordReachedCutoff(indexer.wordInts[whichWord]^xorMask))
    {
        if (whichWord==indexer.wordInts.size()-1)
        {
            count = wordHashLeft->findGenomeLocations(indexer.wordInts[whichWord]^xorMask, indexer.wordInts[whichWord-1], candidates);
        }
        else
        {
            count = wordHashRight->findGenomeLocations(indexer.wordInts[whichWord]^xorMask, indexer.wordInts[whichWord+1], candidates);
        }
    }
    else
    {
        count = wordIndex->hashindices[(indexer.wordInts[whichWord]^xorMask) + 1] - wordIndex->hashindices[indexer.wordInts[whichWord]^xorMask];
        candidates = &wordIndex->wordpositions[wordIndex->hashindices[indexer.wordInts[whichWord]^xorMask]];
    }

#if  0
    // check how many count and where the word matched positions
    if (count > 0)
    {
        std::cerr << "count = " << count << " candidates = ";
        for (int i = 0; i< count; i++)
            std::cerr<< candidates[i] << ",";
        std::cerr<<" at " << __FILE__ << ":" << __LINE__ << std::endl;
    }
#endif

    if (allPositionsMethod)
    {
        return (*allPositionsMethod)(
                   this,
                   indexer,
                   count,
                   candidates,
                   whichWord
               );
    }
    else
    {
        return evalAllCandidatePositions(
                   onePositionMethod,
                   indexer,
                   whichWord,
                   count,
                   candidates
               );
    }

    return false;
}

// a wrapper function
// to call either "allPositionMethod" (a function pointer passed in as parameter)
// or calling "evalAllCandidatePositions" function using parameter "onePositionMethod"
bool MapperBase::evalAllCandidatesForColorSpaceWord(
    evalSinglePositionFunctionType onePositionMethod,
    evalAllPositionsFunctionType allPositionsMethod,
    ReadIndexer &indexer,
    unsigned int whichWord,
    wordInteger_t shiftLocation,
    wordInteger_t mutationIndex)
{

    wordInteger_t mask = 15;
    mask<<= shiftLocation;
    unsigned int adjacentColorCodes = (indexer.wordInts[whichWord] & mask) >>shiftLocation;
    wordInteger_t word = (indexer.wordInts[whichWord] & ~mask) |
                         (colorSpaceSNP[adjacentColorCodes][mutationIndex] << shiftLocation);
    word &= wordIndex-> wordsCountIndexMask; // cap "word" to its maximum legal value

    int count = 0;
    genomeIndex_t *candidates = NULL;

    if (wordIndex->wordReachedCutoff(word))
    {
        if (whichWord==indexer.wordInts.size()-1)
        {
            count = wordHashLeft->findGenomeLocations(word, indexer.wordInts[whichWord-1], candidates);
        }
        else
        {
            count = wordHashRight->findGenomeLocations(word, indexer.wordInts[whichWord+1], candidates);
        }
    }
    else
    {
        count = wordIndex->hashindices[(word) + 1] - wordIndex->hashindices[word];
        candidates = &wordIndex->wordpositions[wordIndex->hashindices[word]];
    }

    if (allPositionsMethod)
    {
        return (*allPositionsMethod)(
                   this,
                   indexer,
                   count,
                   candidates,
                   whichWord
               );
    }
    else
    {
        return evalAllCandidatePositions(
                   onePositionMethod,
                   indexer,
                   whichWord,
                   count,
                   candidates
               );
    }

    return false;
}

inline bool MapperBase::evalAllCandidatePositions(
    evalSinglePositionFunctionType pmethod,
    ReadIndexer &indexer,
    int     whichWord,
    int     candidateCount,
    genomeIndex_t *candidates)
{
    for (int i = 0; i < candidateCount; i ++)
    {
#if defined(LIMIT_CANDIDATES)
        if (++totalCandidateCount>LIMIT_CANDIDATES) return true;
#endif
        // avoid underflow in next subtract:
        if (indexer.wordPositions[whichWord] > candidates[i]) continue;

        // here, candidate position can (rarely) be smaller than the read
        // position ... see Test.cpp for test case
        //
        genomeIndex_t genomeMatchPosition =  candidates[i] - indexer.wordPositions[whichWord];

        // if genomeMatchPosition goes beyond the whole length of the genome
        if ((genomeMatchPosition + indexer.read.size()) > gs->sequenceLength())
            continue;

        // for color space, since we chopped the first two base, we should check the beginning part of the sequence
        if (genomeMatchPosition < 2 && indexer.gs->isColorSpace())
            continue;

        // check if this position has already been evalulated;
        if (indexer.checkedPositions.Find(genomeMatchPosition) != -1)
            continue;

        indexer.checkedPositions.Add(genomeMatchPosition);
        if ((*pmethod)(this, indexer, genomeMatchPosition, whichWord)) return true;
    }
    return false;
}
