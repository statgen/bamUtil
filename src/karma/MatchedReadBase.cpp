#include <cmath>
#include "Debug.h"
#include "../bam/SamFlag.h"

#include "MatchedReadBase.h"
#include "ColorSpace.h"

static inline char baseComplement(const char& c) {
    return GenomeSequence::base2complement[(int) c];
};

#include "MatchedReadBase.h"
#include "ColorSpace.h"
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
    genomeMatchPosition = INVALID_GENOME_INDEX;
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
                    // read bases but not the reference
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
                    // programming error, as the code should not reach here.
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
} //bool MatchedReadBase::getBSSeq2print(std::string&    sequence2print, ...


/*
 *  Explanatin of a position relationship graph
 *  Forward match:
 *  reference index:        0 1 2 3 4 5 6 7 8
 *  reference (BS):         A C T G A C G T A
 *  reference (CS):           1 2 1 2 1 3 1 3
 *  read (CS):            A 0|1 2 1 2               (note: 1 2 1 2 is the part that actually aligned)
 *  read (BS):           (A)A|C T G A               (note: (A) is a primer. it should not be outputed)
 *  read index(CS):       0 1|2 3 4 5
 *  genomeMatchPosition:     |^                     (note: ^ marked the genomeMatchPosition)
 *  E.g. cs_index = 3 => genome_index = referencePosition - 2 + (3) = referencePosition + 1
 *
 *  Backward match:
 *  reference index:        0 1 2 3 4 5 6 7 8
 *  reference (BS):         A C T G A C G T A
 *  reference (CS):           1 2 1 2 1 3 1 3
 *  read (CS):                        1 3 1 3|0 T     (note: 1 3 1 3 is the part that actually aligned)
 *  read (BS - observed.)             T G C A|T(T)    (note: (T) is a primer, should not be outputed)
 *  read (BS - rev. comp.):           A C G T|A(A)    (note: reverse complement in base space, it matches BS from index 4)
 *  read index (CS):                  5 4 3 2|1 0
 *  genomeMatchPosition:              ^               (note: ^ marked the genomeMatchPosition)
 *  E.g. cs_index = 2 => genome_index = referencePosition + (6-1) - (2) = referencePosition + 3
 *
 *  NOTE: 
 *  when aligning the color space reads, there are two indices here, one is the CS read index,
 *  the other is reference index (of the type genomeIndex_t) used in base space and color space sequence.
 *  we assume the base space (BS) read that will be translated frmo CS read will have the same 
 *  length as the CS read, and only in the last step, we will chop off the first leading primer.
 *
 *  in the FORWARD match case, the basef (CS) at poisition i (read index) represents
 *  a trnsiation of the BS read at positino (i-1) and (i), that corresponds to the 
 *  reference index (genomeMatchPosition + i -2) 
 *  an example is the '2' above the '^'. In read index = 2 corresponds to the genomeMatchPosition = 1, 
 *  we can verify that (genomeMatchPosition + 2 - 2) = genomeMatchPosition
 *  We can also verify that the color base = '1' at that position represents a transition of 
 *  the change in BS at index (2-1) = 1 and (2) = 2 , aka the change of 'A' to 'C' is color '1'. 
 *
 *  in the BACKWARD match case, the base (CS) at poisition i (in read index) also represents
 *  a trnsiation of the BS read at positino (i-1) and (i), that corresponds to the 
 *  reference index (genomeMatchPosition + length -1 ) - (i-1) and (genomeMatchPosition + length -1 ) - (i)
 *  an example is the '5' above the '^'. In read index = 5 corresponds to the genomeMatchPosition = 5
 *  we can verify that (genomeMatchPosition + 6 - 1) - (5) = genomeMatchPosition
 *  we can also verify that the color base = '1' at that position represents a transition of 
 *  the change in BS at index (5-1) = 4 and (5), aka. the change of 'G' to 'T' is color '1'.
 *  
 *  Another issue of the strand direction is that the translated BS read may not corresponds to 
 *  the BS genome sequence in the same way. For example, in the forward strand match, as in the first graph,
 *  the translated BS read matches the BS genome seuqence. 'CTGA' of the read are in the same column of 
 *  the BS genome sequence. But in the forward strand match, as in the second graph,
 *  the translated BS read only matches BS genome sequence with a OFFSET of 1. 
 *  'ACGT' of the complement of the read begins at reference index = 5 but it is aligned to the 
 *  reference index = 4. This minor difference bring many pitfalls when we translating the CS read to BS read.
 */
void MatchedReadBase::translateCS2BSByCigarOperation(const Cigar::CigarOperator& cigar,
                                                     std::string& cs_read,
                                                     std::string& cs_qual,
                                                     uint32_t& readPosition,
                                                     GenomeSequence* gs,
                                                     GenomeSequence* csgs,
                                                     genomeIndex_t& referencePosition,
                                                     bool isForward,
                                                     std::string& sequence2print,
                                                     std::string& quality2print)
{
    assert(readPosition >= 2);

    uint32_t count = cigar.count;
    if (count == 0) return;

    switch (cigar.operation)
    {
    case CigarRoller::match:
    case CigarRoller::mismatch:
        translateCigarMatchSequence(count,
                                    cs_read, cs_qual, readPosition,
                                    gs, csgs, referencePosition,
                                    isForward,
                                    sequence2print, quality2print );

        readPosition += count;
        if (isForward)
            referencePosition += count;
        else
            referencePosition -= count;
        break;
    case CigarRoller::insert:
        // for an insert, we just translate the bases from last known translated base.
        for (uint32_t i = readPosition ; i < readPosition + count; i++) {
            BaseAndColorToBase( sequence2print[i - 1],
                                cs_read[i],
                                sequence2print[i]);
            quality2print[i] = cs_qual[i];
        }
        // the bases in the read are extra, so we skip over
        // them but not the reference
        readPosition += count;
        break;
    case CigarRoller::del:
        // for delete, there is nothing to print, but we skip reference bases
        if (isForward)
            referencePosition += count;
        else
            referencePosition -= count;
        break;
    case CigarRoller::softClip:
        // for a soft clip, we will translate from last known bases
        // (in other word, we will not the information from the sequence).
        for (uint32_t i = readPosition ; i < readPosition + count; i++) {
            BaseAndColorToBase( sequence2print[i - 1],
                                cs_read[i],
                                sequence2print[i]);
            quality2print[i] = cs_qual[i];
        }
        readPosition += count;
        if (isForward)
            referencePosition += count;
        else
            referencePosition -= count;
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
} // void MatchedReadBase::translateCS2BSByCigarOperation(const Cigar::CigarOperator& cigar, ...)


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

    // if showReference is true:  match will be readBase, mismatch will be lower case of readBase
    // if showReference is false: match will be  "="    , mismatch will be readBase
    if (qualityIsValid())
    {
        assert(genomeMatchPosition != INVALID_GENOME_INDEX );

        int length = cs_read_fragment.size();
        sequence2print.resize(length);
        quality2print.resize(length);
        sequence2print[0] = cs_read_fragment[0];
        quality2print[0] = cs_data_quality[0];
        fill(sequence2print.begin()+1, sequence2print.begin()+length, 'N');
        fill(quality2print.begin()+1, quality2print.begin()+length, '!');

        // we try to synchronize readPosition and referencePosition
        // however, remember that we will have to manually recover the base after the primer
        uint32_t readPosition = 2;
        genomeIndex_t referencePosition = genomeMatchPosition;
        if (!indexer->isForward) 
        {
            referencePosition = (genomeMatchPosition + length - 1) - (2); 
        }

        // need to rewrite calibrateSequence() function as we need to support cigar
        // the work flow is:
        // follow the CIGAR, translate CS to BS sequence by one CIGAR operator at a step
        // when CIGAR operator == M, we may need to call fixBaseRange() function
        // when CIGAR operator == S, we just skip it
        // when CIGAR operator == I/D, we just change the position
        // due to the complexity of the translation process, we will handle forward
        // and backward separately

        if (indexer->isForward) // forward match
        {
            // recover the first base
            BaseAndColorToBase(cs_read_fragment[0], cs_read_fragment[1], sequence2print[1]);
            quality2print[1] = cs_data_quality[1];
            if (sequence2print[1] == (*gs)[genomeMatchPosition-1])
            {
                CigarRoller temp("1M");
                temp+= cigarRoller;
                cigarRoller = temp;
            } else {
                // we can set add '1S' to the leftmost of the cigar,
                // however, the next step, translateCS2BSByCigarOperation() will be affetced.
                // for simplicity we will also append '1M' here for now.
                CigarRoller temp("1M");
                temp+= cigarRoller;
                cigarRoller = temp;
            }

            // recover the rest
            for (int cigarIndex = 0; cigarIndex < cigarRoller.size(); cigarIndex++)
            {

                // using cigarRoller[cigarIndex].operation, cigarRoller[cigarIndex].count
                // to translate a piece of color space code
                // the translated piece of bases will be concatnated to
                // sequence2print and quality2print
                translateCS2BSByCigarOperation(cigarRoller[cigarIndex],
                                               cs_read_fragment, cs_data_quality, readPosition,
                                               gs, csgs, referencePosition,
                                               indexer->isForward,
                                               sequence2print, quality2print);
            }
        } else{ // backward match
            // recover the first base
            BaseAndColorToBase(cs_read_fragment[0], cs_read_fragment[1], sequence2print[1]);
            quality2print[1] = cs_data_quality[1];

            // recover the rest
            for (int cigarIndex = 0; cigarIndex < cigarRoller.size(); cigarIndex++)
            {

                // using cigarRoller[cigarIndex].operation, cigarRoller[cigarIndex].count
                // to translate a piece of color space code
                // the translated piece of bases will be concatnated to
                // sequence2print and quality2print
                translateCS2BSByCigarOperation(cigarRoller[cigarIndex],
                                               cs_read_fragment, cs_data_quality, readPosition,
                                               gs, csgs, referencePosition,
                                               indexer->isForward,
                                               sequence2print, quality2print);
            }

            // handle CIGAR for the base after the primer
            if (sequence2print[1] == baseComplement( (*gs)[(genomeMatchPosition+length-1)-(1)] ) )
            {
                cigarRoller.Add(CigarRoller::match, 1);
            } else {
                //cigarRoller.Add(CigarRoller::match, 1);
                cigarRoller.Add(CigarRoller::softClip, 1);
            }
        } // end if (indexer->isForward)

        // chop the first base
        sequence2print.erase(0, 1);
        quality2print.erase(0, 1);

        if (!indexer->isForward) {
            std::reverse(sequence2print.begin(), sequence2print.end());
            //sequence2print = getComplement(sequence2print);
            std::reverse(quality2print.begin(), quality2print.end());
        }

        // translate matched bases to = or lower case.
        if (indexer->isForward)
            // -1 is because we want to begin comparison from the base before the genomeMatchPosition
            // e.g. (A) T A C A C       (A) is the primer, here sequence2print does not contain A
            //            ^             the genomeMatchPosition
            // we want to start comparison from the 'T'
            markUnmatchedBases(cigarRoller, gs, genomeMatchPosition-1, showReferenceBases, sequence2print);
        else 
            // for backward strand, things are simpler, since the genomeMatchPosition is the leftmost part 
            // e.g. C A C A T (A)       (A) is the primer, here sequence2print does not contain A
            //      ^                   the genomeMatchPosition
            markUnmatchedBases(cigarRoller, gs, genomeMatchPosition, showReferenceBases, sequence2print);
    }
    else
    {
        // if the quality is invalid, trying to compare against reference is useless,
        // so at least just print the read here...
        sequence2print = convertCSRead(cs_read_fragment);
        quality2print = convertCSQuality(cs_data_quality);
    }
    return true;
} //bool MatchedReadBase::getCSSeqAndQual2print(std::string&    sequence2print,...)


/*
 * mark unmatched bases by showReferenceBases
 *   if (showReferenceBases == true),  match will be marked as A, mismatch will be marked as a
 *   if (showReferenceBases == false), match will be marked as =, mismatch will be marked as A
 *   for indels, all bases will be shown like mismatches.
 * NOTE: sequence2print does NOT contain the leading primer
 * NOTE: we will work on FORWARD strand in base space.
 * NOTE: we will NOT check the length given by CigarRoller equal to the length of sequence2print
 */
void MatchedReadBase::markUnmatchedBases(const CigarRoller& cigarRoller, 
                                         GenomeSequence* gs, genomeIndex_t genomePos, 
                                         bool showReferenceBases,
                                         std::string& sequence2print)
{
    assert(genomePos != INVALID_GENOME_INDEX);

    uint32_t readPos = 0;
    for (int cigarIndex = 0; cigarIndex < cigarRoller.size(); ++cigarIndex) {
        const Cigar::CigarOperator& cigar = cigarRoller[cigarIndex];
        if (cigar.count == 0) continue;
        switch (cigar.operation) {
        case CigarRoller::match:
        case CigarRoller::mismatch:
            for (uint32_t i = 0; i < cigar.count; i++) {
                if (showReferenceBases == true) {
                    if (sequence2print[readPos] == (*gs)[genomePos]) {
                        // we don't need to change anything (it's already like A)
                    } else {
                        sequence2print[readPos] = tolower(sequence2print[readPos]);
                    }
                } else { // do not show reference bases 
                    if (sequence2print[readPos] == (*gs)[genomePos]) {
                        sequence2print[readPos] = '=';
                    } else{
                        // we don't need to change anything. (it's already like A)
                    }
                }
                ++ readPos; 
                ++ genomePos;
            }
            break;
        case CigarRoller::del:
            genomePos += cigar.count;
            break;
        case CigarRoller::insert:
            for (uint i = 0; i < cigar.count; i++) {
                if (showReferenceBases == true) 
                {
                    sequence2print[readPos] = tolower(sequence2print[readPos]);
                } else {
                    // do not need to consider this case, as the bases are already upper cases
                }
                readPos ++;
            }
            break;
        case CigarRoller::softClip:
            for (uint i = 0; i < cigar.count; i++) {
                if (showReferenceBases == true) {
                    if (sequence2print[readPos] == (*gs)[genomePos]) {
                        // we don't need to change anything (it's already like A)
                    } else {
                        sequence2print[readPos] = tolower(sequence2print[readPos]);
                    }
                } else { // do not show reference bases 
                    if (sequence2print[readPos] == (*gs)[genomePos]) {
                        sequence2print[readPos] = '=';
                    } else{
                        // we don't need to change anything. (it's already like A)
                    }
                }
                readPos ++;
                genomePos ++;
            }
            break;
        case CigarRoller::hardClip:
        case CigarRoller::pad:
            break;
        case CigarRoller::none:
        default:
            // unsupported Cigar operatio goes here.
            assert(false);
            break;
        }
    }
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
// QNAME    [^ \t\n\r]+             Query pair NAME if paired; or Query NAME if unpaired 2
// FLAG     [0-9]+ [0,216-1]        bitwise FLAG (Section 2.2.2)
// RNAME    [^ \t\n\r@=]+           Reference sequence NAME 3
// POS      [0-9]+    [0,229-1]     1-based leftmost POSition/coordinate of the clipped sequence
// MAPQ     [0-9]+    [0,28-1]      MAPping Quality (phred scaled prob. that the alignment is wrong) 4
// CIGAR    ([0-9]+[MIDNSHP])+|\*   extended CIGAR string
// MRNM     [^ \t\n\r@]+            Mate Reference sequence NaMe; “=” if the same as <RNAME> 3
// MPOS     [0-9]+    [0,229-1]     1-based leftmost Mate POSition of the clipped sequence
// ISIZE    -?[0-9]+    [-229,229]  inferred Insert SIZE 5
// SEQ      [acgtnACGTN.=]+|\*      query SEQuence; “=” for a match to the reference; n/N/. 
//                                  for ambiguity; cases are not maintained 6,7
// QUAL     [!-~]+|\*               query QUALity; ASCII-33 gives the Phred base quality 6,7
// TAG      [A-Z][A-Z0-9]           TAG
// VTYPE    [AifZH]                 Value TYPE
// VALUE    [^\t\n\r]+              match <VTYPE> (space allowed)
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
void MatchedReadBase::print(std::ostream      &file,
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

    const char *chromosomeName = "*";
    int chromosomePosition = 0;
    // for paired reads, set the MRNM field.
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
    if (indexer->phredQuality!="*") assert(sequence2print.size() == indexer->phredQuality.size());

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
        insertSize = (int64_t) mate->genomeMatchPosition - (int64_t) genomeMatchPosition ;
        if (insertSize >= 0) // mate coordinate is bigger
        {
            // add mate length (no matter mate's strand)
            insertSize += mate->indexer->read.size();
        }
        else // mate coordinate is smaller, insertSize < 0
        {
            // minus this length (no matter mate's strand)
            insertSize -= this->indexer->read.size();
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

/**
 * Fix a range of mismatches in color space
 * and assign corresponding qualities for those mismatched bases
 * (that could happend when translating color space reads to base space)
 *
 * for 1 consecutive mismatch, we will translate that mismatched base to most likely one
 * for 2 consecutive mismatches, we may be able to translate it as a SNP or 
 * pick alternative bases according to qualities
 * for more consecutive mismatches, we will have to skip here. However, we expect this is the very rare case.
 *
 * One example of why we are doing this is:
 * start = 1, end = 2   (inclusive for start and end, forward strand
 *                       they are read index (CS) as shown in the following diagram)
 * position:   0123....
 * color read:  00  (position 1, 2)
 * ref genome: ATCC
 * ref genome: ?320
 * translate1: AAAA (wrong, as it makes long mismatches in base space)
 * translate2: AACC (correct, it makes use of the information at poisition 0 and 3 of the ref genome)
 * 2 mismatches, at position 1 and 2, translate2 is preferred over translate1.
 * Detailed algorithm is shown in the function body.
 *
 * NOTE: referencePosition does not equal to genomeMatchPosition, instead, it tries to match the readPos
 */
void MatchedReadBase::fixBaseRange(uint32_t start, uint32_t end, /// inclusive boundaries of the BS read
                                   std::string& read_fragment, std::string& data_quality,
                                   const std::string& cs_read_fragment, const std::string& cs_data_quality,
                                   GenomeSequence* gs, GenomeSequence* csgs,
                                   genomeIndex_t referencePosition,
                                   bool isForwardStrand)
{
    if (referencePosition == INVALID_GENOME_INDEX) return;
    // we can only fix mismatches after start =2
    assert(start>=2);

    // we will not fix longer (>2) color space mismatches
    // as single color space error rate is 1-2%, at least 3 consecutive error has rate 1^(-6), which is rare
    int nMismatch = end-start+1;
    if (nMismatch>2)
        return;

    // fix the reads end part, as it is usually easier.
    if (end == read_fragment.size() - 1) {
        if (isForwardStrand) {
            for (uint32_t i = start ; i <= end; i++) {
                // TODO(zhanxw): will handle cs_read_fragment[i] == '5' case more gracefully.
                BaseAndColorToBase(read_fragment[i-1], cs_read_fragment[i], read_fragment[i]);
                data_quality[i] = cs_data_quality[i];
            }
        } else { // backward strand
            for (uint32_t i = start ; i <= end; i++) {
                // TODO(zhanxw): will handle cs_read_fragment[i] == '5' case more gracefully.
                BaseAndColorToBase(read_fragment[i-1], cs_read_fragment[i], read_fragment[i]);
                data_quality[i] = cs_data_quality[i];
            }
        }
        return;
    }

    // for 1 mismatched BS base, there are TWO color bases use its information
    // thus we can infer the mismatched base call in two ways.
    // for 2 mismatched BS bases, there are THREE color bases involved.
    // if the 2 mismatched BS bases the results of SNP
    // we will pick the highest TWO color bases and use them to infer the BS bases.
    char choice1, choice2;
    char qual1, qual2, qual3;

    if (isForwardStrand)
    {
        switch (nMismatch)
        {
        case 1:
            //  Method 1:
            //    ___________
            //   |     |     |
            //   |  X --> ?  |              (BS genome seqeunce)
            //   |_____|_____| 
            //          start               (CS read, we will use cs_read[start] and X to infer ? , 
            //           end                 since the ? in BS genome sequence and 
            //    ___________                translated BS read are the same)
            //   |     |     |
            //   |     |  ?  |              (to be translated BS read)
            //   |_____|_____| 
            //
            //  Method 2:
            //    ___________
            //   |     |     |
            //   |  ? <-- X  |              (BS genome seqeunce)
            //   |_____|_____| 
            //    start                     (CS read, we will use cs_read[start] and X to infer ? , 
            //     end                       since the ? in BS genome sequence and 
            //    ___________                translated BS read are the same)
            //   |     |     |
            //   |  ?  |     |              (to be translated BS read)
            //   |_____|_____| 
            //
            // here choice1 is using Method 1
            BaseAndColorToBase((*gs)[referencePosition-2 + start-1], cs_read_fragment[start], choice1);
            // here choice1 is using Method 2
            BaseAndColorToBase((*gs)[referencePosition-2 + (end+1)], cs_read_fragment[end+1], choice2);
            qual1 = data_quality[(start)];
            qual2 = data_quality[(end+1)];
            if (choice1 != choice2)
            {
                // pick base (BS) by higher quality, reduce its quality
                // by empirically formlar: high_base_quality - low_base_quality
                if (qual1 > qual2)
                {
                    read_fragment[(start)] = choice1;
                    data_quality[(start)] = qual1 - qual2 + 33;
                }
                else
                {
                    read_fragment[(start)] = choice2;
                    data_quality[(start)] = qual2 - qual1 + 33;
                }
            }
            else
            {
                read_fragment[(start)] = choice1 ;
                // set to lower quality
                if (qual1 > qual2)
                    data_quality[(start)] =qual2;
                else
                    data_quality[(start)] =qual1;
            }
            break;
        case 2:
            //    ___________ ___________ 
            //   |     |     |     |     |
            //   |  X --> ?  |  * <-- Y  |      (BS genome seqeunce)
            //   |_____|_____|_____|_____| 
            //                                  (CS read, we will use cs_read[start] and X to infer ?, and
            //           start end               use cs_read[end+1] and Y to infer *, then check if 
            //    ___________ ___________        the transition of ? to * equal to cs_read[end]
            //   |     |     |     |     |
            //   |     |  ?  |  *  |     |      (to be translated BS read)
            //   |_____|_____|_____|_____| 
            //
            char color;
            BaseAndColorToBase((*gs)[(referencePosition-2)+(start-1)], cs_read_fragment[start], choice1);
            BaseAndColorToBase((*gs)[(referencePosition-2)+(end+1)], cs_read_fragment[end+1], choice2);
            BaseAndBaseToColor(choice1, choice2, color);
            if (color != cs_read_fragment[end])   // Not a SNP, so it likely there's one wrong color base.
            {
                // we pick the two highest qualities surrounding 2 mismatch position
                // the use these 2 colors to call bases
                qual1 = data_quality[(start)];
                qual2 = data_quality[(end)];
                qual3 = data_quality[(end+1)];
                if (qual1 >= qual3 && qual2 >= qual3)   //qual3 is the smallest
                {
                    read_fragment[(start)] = choice1;
                    BaseAndColorToBase(choice1, cs_read_fragment[end], read_fragment[(end)]);
                    data_quality[start] = qual1 - qual3 +33;
                    data_quality[end] = qual2 - qual3 +33;
                }
                else if (qual1 >= qual2 && qual3 >= qual2)   //qual2 is the smallest
                {
                    read_fragment[(start)] = choice1;
                    read_fragment[(end)] = choice2;
                    data_quality[start] = qual1 - qual2 +33;
                    data_quality[end] = qual3 - qual2 +33;
                }
                else   //qual1 is the smallest
                {
                    read_fragment[(end)] = choice2;
                    BaseAndColorToBase(choice2, cs_read_fragment[end], read_fragment[(start)]);
                    data_quality[start] = qual2 - qual1 +33;
                    data_quality[end] = qual3 - qual1 +33;
                }
            }
            else    // 2 consecutive msimatch color codes represents a SNP
            {
                read_fragment[(start)] = choice1;
                read_fragment[(end)] = choice2;
                data_quality[(start)] = cs_data_quality[start];
                data_quality[(end)] = cs_data_quality[end];
            }
            break;
        default:
            return;
        }
    }
    else   // backwardStrand
    {
        char base1 ='N', base2='N';
        switch (nMismatch)
        {
        case 1:
            //  Method 1:
            //    ___________
            //   |     |     |
            //   |  ? <-- X  |                  (BS genome seqeunce)
            //   |_____|_____| 
            //          start                   (CS read, we will use cs_read[start] and X to infer ? , 
            //           end                     since the ? in BS genome sequence and 
            //          ___________              translated BS read are the same)
            //         |     |     |
            //         |  ?  |     |            (to be translated BS read)
            //         |_____|_____| 
            //
            //  Method 2:
            //    ___________
            //   |     |     |
            //   |  X --> ?  |                  (BS genome seqeunce)
            //   |_____|_____| 
            //                start             (CS read, we will use cs_read[end+1] and X to infer ? , 
            //                 end               since the ? in BS genome sequence and 
            //          ___________              translated BS read are the same)
            //         |     |     |
            //         |     |  ?  |            (to be translated BS read)
            //         |_____|_____| 
            //
            // here choice1 is using Method 1
            base1 = baseComplement( (*gs)[(referencePosition)-(start)] ); 
            BaseAndColorToBase(base1, cs_read_fragment[start], choice1);
            // here choice1 is using Method 2
            base2 = baseComplement( (*gs)[(referencePosition)-(end+2)] );
            BaseAndColorToBase(base2, cs_read_fragment[end+1], choice2);
            qual1 = data_quality[(start)];
            qual2 = data_quality[(end+1)];
            if (choice1 != choice2)
            {
                if (qual1 > qual2)
                {
                    read_fragment[(start)] = choice1;
                    data_quality[(start)] = qual1 - qual2 + 33;
                }
                else
                {
                    read_fragment[(start)] = choice2;
                    data_quality[(start)] = qual2 - qual1 + 33;
                }
            }
            else
            {
                read_fragment[(start)] = choice1 ;
                if (qual1 > qual2)
                    data_quality[(start)] =qual2;
                else
                    data_quality[(start)] =qual1;
            }
            break;
        case 2:
            //    ___________ ___________ 
            //   |     |     |     |     |
            //   |  X --> ?  |  * <-- Y  |                  (BS genome seqeunce)
            //   |_____|_____|_____|_____| 
            //                                              (CS read, we will use cs_read[start] and Y to infer * , 
            //                 end  start                    and use cs_read[end+1] and X to infer ?, and check if the
            //                                               transition of * to ? equal to cs[end] )
            //          ___________ ___________              
            //         |     |     |     |     |             
            //         |     |  ?  |  *  |     |            (to be translated BS read)
            //         |_____|_____|_____|_____| 
            //
            char color;
            base1 = baseComplement( (*gs)[(referencePosition)-(start)] );
            base2 = baseComplement( (*gs)[(referencePosition)-(end+2)] );

            BaseAndColorToBase(base1, cs_read_fragment[start], choice1);
            BaseAndColorToBase(base2, cs_read_fragment[end+1], choice2);
            BaseAndBaseToColor(choice1, choice2, color);
            if (color != cs_read_fragment[end])
            {
                // we pick the two highest qualities surrounding 2 mismatch position
                // the use these 2 colors to call bases
                qual1 = data_quality[(start)];
                qual2 = data_quality[(end)];
                qual3 = data_quality[(end+1)];
                if (qual1 >= qual3 && qual2 >= qual3)   //qual3 is the smallest
                {
                    read_fragment[(start)] = choice1;
                    BaseAndColorToBase(choice1, cs_read_fragment[end], read_fragment[(end)]);
                    data_quality[start] = qual1 - qual3 +33;
                    data_quality[end] = qual2 - qual3 +33;
                }
                else if (qual1 >= qual2 && qual3 >= qual2)   //qual2 is the smallest
                {
                    read_fragment[(start)] = choice1;
                    read_fragment[(end)] = choice2;
                    data_quality[start] = qual1 - qual2 +33;
                    data_quality[end] = qual3 - qual2 +33;
                }
                else   //qual1 is the smallest
                {
                    read_fragment[(end)] = choice2;
                    BaseAndColorToBase(choice2, cs_read_fragment[end], read_fragment[(start)]);
                    data_quality[start] = qual2 - qual1 +33;
                    data_quality[end] = qual3 - qual1 +33;
                }
            }
            else    // 2 consecutive msimatch color codes represents a SNP
            {
                read_fragment[(start)] = choice1;
                read_fragment[(end)] = choice2;
                data_quality [(start)] = cs_data_quality[start];
                data_quality [(end)] = cs_data_quality[end];
            }
            break;
        default:
            return;
        }
    }
} // void MatchedReadBase::fixBaseRange(int start, int end, /// inclusive boundaries...

void MatchedReadBase::translateCigarMatchSequence(uint32_t count,
                                                  std::string& cs_read, std::string& cs_qual, int readPosition,
                                                  GenomeSequence* gs, GenomeSequence* csgs, genomeIndex_t referencePosition,
                                                  bool isForward,
                                                  std::string& sequence2print, std::string& quality2print)
{
    bool lastBaseMatched = true;
    int start = 0;
    int end = 0;
    if (isForward) {
        for (uint32_t i = 0; i < count; i++)
        {
            int readPos = readPosition +i ;
            int genomePos = referencePosition + i ;
            if (cs_read[readPos] == (*csgs)[genomePos] ) { // this bases match
                if (lastBaseMatched == true) {
                    // although we can use CS translation to assign
                    // this base, just for simplicity, we use gs (BS) directly.
                    sequence2print[readPos] = (*gs)[genomePos];
                    quality2print[readPos] = cs_qual[readPos];
                }
                if (lastBaseMatched == false) {
                    lastBaseMatched = true;
                    end = readPos-1;
                    if (end != 0) {
                        // fix [start, end] bases
                        fixBaseRange(start, end, sequence2print, quality2print,
                                     cs_read, cs_qual, gs, csgs, referencePosition, isForward);
                    }
                    // fix (end+1) base, aka (readPos) base
                    BaseAndColorToBase((*gs)[genomePos-1], cs_read[readPos], sequence2print[readPos]);
                }
                
            } else { // does not match
                if (lastBaseMatched == true) {
                    lastBaseMatched = false;
                    start = readPos;
                }  else {
                    end = readPos;
                    fixBaseRange(start, end, sequence2print, quality2print,
                                 cs_read, cs_qual, gs, csgs, referencePosition, isForward);
                }
            }
        }
        if (lastBaseMatched == false) {
            end = readPosition + count - 1;
            fixBaseRange(start, end, sequence2print, quality2print,
                         cs_read, cs_qual, gs, csgs, referencePosition, isForward);
        }
    } else { // backward match
        for (uint32_t i = 0; i < count; i++)
        {
            int readPos = readPosition +i ;
            int genomePos = referencePosition - i ;
            if (cs_read[readPos] == (*csgs)[genomePos] ) { // this bases match
                if (lastBaseMatched == true) {
                    // although we can use CS translation to assign
                    // this base, just for simplicity, we use gs (BS) directly.
                    sequence2print[readPos] = (*gs)[genomePos];
                    quality2print[readPos] = cs_qual[readPos];
                }
                if (lastBaseMatched == false) {
                    lastBaseMatched = true;
                    end = readPos-1;
                    if (end != 0) {
                        fixBaseRange(start, end, sequence2print, quality2print,
                                     cs_read, cs_qual, gs, csgs, referencePosition, isForward);
                    }
                    char base = baseComplement( (*gs)[genomePos] );
                    // (*gs)[genomePos]  and cs_read[readPos] determines (*gs)[genomePos-1],
                    // which is sequence2print[readPos]
                    BaseAndColorToBase(base, cs_read[readPos], sequence2print[readPos]);
                }

            } else { // does not match
                if (lastBaseMatched == true) {
                    lastBaseMatched = false;
                    start = readPos;
                }  else {
                    end = readPos;
                    fixBaseRange(start, end, sequence2print, quality2print,
                                 cs_read, cs_qual, gs, csgs, referencePosition, isForward);
                }
            }
        }
        if (lastBaseMatched == false) {
            end = readPosition + count - 1;
            fixBaseRange(start, end, sequence2print, quality2print,
                         cs_read, cs_qual, gs, csgs, referencePosition, isForward);
        }
    }
} //void MatchedReadBase::translateCigarMatchSequence(uint32_t count,...)

void MatchedReadBase::printColorSpace(std::ostream &file,
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
        cigarString = cigarRoller.getString();
    }
    else
    {
        cigarString = defaultCigarString;
    }

    //
    // final sanity check:
    //  read length == qual length == cigar string M count
    if (indexer->phredQuality!="*") assert(sequence2print.size() == quality2print.size());
    if (strcmp(cigarString, "*")!=0) assert(sequence2print.size() == (unsigned int) cigarRoller.getExpectedQueryBaseCount());

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

void MatchedReadBase::printOptionalTags(std::ostream &file, 
                                        bool isPairAligned, 
                                        const std::string &sampleGroupID, 
                                        const std::string &alignmentPathTag)
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
    std::cout << (isForward() ? 'F' : 'R')
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

#if 0 

//////////////////////////////////////////////////////////////////////
// OBSOLETE COMMENTS
//////////////////////////////////////////////////////////////////////
/***
 * When translate color space reads to base space, it is possible we have long strands of mismatch
 * therefore, we will try to mark the start and end position for the consecutive mismatches,
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
        read_fragment[length-1-1] = baseComplement( (read_fragment[length-1-1]) );
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

// when we access CS read at position (0-based) i, then (shown above, read in the column)
// we are accessing BS read at position (0-based) i-1, and
// the genome position is (genomeMatchPosition+(i)-2)
// e.g. read index(CS) = 2 =>
//      read index(BS) = 2 - 1 = 1 =>
//      genome position = genomeMatchPosition + 2 -2 = genomeMatchPosition
// #define FWD_BSREAD_POS(i) ((i)-1)
// #define FWD_GENOME_POS(i) (genomeMatchPosition+(i)-2)

// when we access CS read at position (0-based) i, then (shown above, read in the column)
// we are accessing BS read at position ((length-1)-(i))
// the genome position is (genomeMatchPosition+(length-1)-(i))
// e.g. read index(CS) = 2 =>
//      read index(BS) = ((6-1)-2) = 3 =>
//      genome position = genomeMatchPosition + (6-1)-2 = genomeMatchPosition + 3
// #define BWD_BSREAD_POS(i) ((length-1)-(i))
// #define BWD_GENOME_POS(i) (genomeMatchPosition+(length-1)-(i))

#define BACKWARD_CS2BS_REF_OFFSET (-1)

#endif 
