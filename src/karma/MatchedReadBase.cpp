#include <cmath>
#include "debug.h"
#include "../bam/SamFlag.h"

#include "MatchedReadBase.h"

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
