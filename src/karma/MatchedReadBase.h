#ifndef _MATCHEDREADBASE_H_
#define _MATCHEDREADBASE_H_

#include "GenomeSequence.h"
#include "ColorSpace.h"
#include "ReadIndexer.h"
#include "CigarRoller.h"
#include <algorithm>

//
// When we get a match, mutated or not, keep records here
//
class MatchedReadBase
{
  public:
    MatchedReadBase()
    {
        constructorClear();
    }
    virtual ~MatchedReadBase();
    void constructorClear();

    static const int UNSET_QUALITY = -1;         // default value
    static const int INVALID_DATA = -2;          // no used yet
    static const int EARLYSTOP_DUPLICATES = -3;  // early stop due to reaching max number of best matches, with all bases matched
    static const int EARLYSTOP_QUALITY = -4;     // early stop due to reaching max posterior quality
    static const int REPEAT_QUALITY = -5;        // not used yet

    // ---------------------------------
    //  Member variables
    // ---------------------------------
    //
    // whichWord is the index word
    // that was used to obtain this match.  When doing Smith Waterman
    // matching, it tells us where to anchor the forward and
    // reverse SW searches.
    //
    // When gappedAlignment is true, we got this match via using Smith
    // Waterman, which means the ::print method needs to correctly
    // print the corresponding cigar string.
    //
    int     whichWord;
    bool    gappedAlignment;    // Smith Waterman was used to align this result

    genomeIndex_t genomeMatchPosition;
    int mismatchCount;
    int quality;
    double cumulativePosteriorProbabilities;
    int numBest;
    int numMatchContributors;

    ReadIndexer *indexer;

    // if "other"  is given, include its quality as
    // a contributor
    double getQualityScore(double sumQ, double posteriorProbability);
    virtual double getQualityScore();
    inline bool qualityIsValid() const
    {
        return quality >= 0;
    }

    inline bool qualityIsFailed() const
    {
        return quality < UNSET_QUALITY;
    }

    inline bool qualityIsUnset() const
    {
        return quality == UNSET_QUALITY;
    }
    bool isForward() const
    {
        return indexer?indexer->isForward:true;
    }

    virtual void updateMatch(MatchedReadBase &);


    // ---------------------------------
    //  SAM output functions:
    // ---------------------------------
    static void printHeader(std::ostream &);

    void calculateNAMEandPOS(const char* &chromosomeName, int& chromosomePosition,
                             const char* &mateChromosomeName, int& mateChromosomePosition,
                             const MatchedReadBase* mate);
    bool getBSSeq2print(std::string&    sequence2print,
                        CigarRoller&    cigarRoller,
                        bool            showReferenceBases = true);
    bool getCSSeqAndQual2print(std::string& sequence2print,
                               std::string& quality2print,
                               std::string&    cs_read_fragment,
                               std::string&    cs_data_quality,
                               GenomeSequence* gs,
                               GenomeSequence* csgs,
                               CigarRoller& cigarRoller,
                               bool         showReferenceBases = true);
    void translateCS2BSByCigarOperation(const Cigar::CigarOperator& cigar, 
                                        std::string& cs_read, std::string& cs_qual, 
                                        uint32_t& readPos, 
                                        GenomeSequence* gs, GenomeSequence* csgs, 
                                        genomeIndex_t& referencePos,
                                        bool showReferenceBases,
                                        std::string& sequence2print, std::string& quality2print);
    void markUnmatchedBases(const CigarRoller& cigarRoller, 
                            GenomeSequence* gs, genomeIndex_t genomePos, 
                            bool showReferenceBases,
                            std::string& sequence2print);
    ///
    /// Write a single read to the stream.  It may
    /// have been sequenced as a pair, in which case
    /// mate points to the mate matchedReadBase.
    ///
    /// outputStream -> ostream that we send output to.
    ///
    /// mate -> if sequenced as a mate, this is the mate.
    /// Having mate!=NULL does not guarantee the mate was
    /// pair aligned (see isProperAligned).
    ///
    /// fragmentTag -> the tag for this read.  Since
    /// it is not included in the indexer class, we have
    /// to include it here.
    ///
    /// cigarRoller -> Vector of CIGAR operations indicating
    /// insert/delete/match/mismatch for the aligned read.
    /// This is populated by either the Smith Waterman aligner
    /// or by some other local alignment functions.
    ///
    /// isProperAligned -> indicates when a mate exists
    /// and is within the max insert distance specified
    /// by the command line arguments.  This affects the
    /// SAM flag bit 0x0002 and ISIZE fields.
    ///
    void print(
        std::ostream &outputStream,
        MatchedReadBase *mate,
        std::string &fragmentTag,
        bool showReferenceBases,
        CigarRoller  &cigarRoller,
        bool    isProperAligned,
        uint16_t    samMateFlag,
        const std::string &sampleGroupID,
        const std::string &alignmentPathTag
        );

    ///
    /// It is a similar procedure to print(), internally, we translate
    /// color space reads to base space, and color space qualities to base space qualities
    ///
    void printColorSpace(
        std::ostream &file,
        GenomeSequence* gs,
        GenomeSequence* csgs,
        MatchedReadBase *mate,
        std::string &cs_read_fragment,
        std::string &cs_data_quality,
        std::string &fragmentTag,
        bool showReferenceBases,
        CigarRoller  &cigarRoller,
        bool    isProperAligned,
        uint16_t    samMateFlag,
        const std::string &sampleGroupID,
        const std::string &alignmentPathTag
        );

    virtual void printOptionalTags(
        std::ostream &, 
        bool isProperAligned,
        const std::string &readGroupID,
        const std::string &alignmentPathTag);


    // MapperBase has debugPrint also, but it does more validation.  This
    // method just prints the data to cout:
    void debugPrint();

    // Fix a range of mismatches in color space
    void fixBaseRange(uint32_t start, uint32_t end, /// inclusive boundaries
                      std::string& read_fragment, std::string& data_quality,
                      const std::string& cs_read_fragment, const std::string& cs_data_quality,
                      GenomeSequence* gs, GenomeSequence* csgs,
                      genomeIndex_t genomeMatchPosition,
                      bool isForwardStrand);

    //
    void translateCigarMatchSequence(uint32_t count, 
                                     std::string& cs_read, std::string& cs_qual, int readPosition,
                                     GenomeSequence* gs, GenomeSequence* csgs, genomeIndex_t referencePosition,
                                     bool isForward, 
                                     std::string& sequence2print, std::string& quality2print) ;


    //////////////////////////////////////////////////
    // Obselete code goes here
    //////////////////////////////////////////////////
#if 0
    void calibrateSequence(std::string& read_fragment, std::string& data_quality,
                           const std::string& cs_read_fragment, const std::string& cs_data_quality,
                           GenomeSequence* gs, GenomeSequence* csgs,
                           genomeIndex_t genomeMatchPosition,
                           bool isForwardStrand);
#endif

#if 0
    // not ever used funcitons
    virtual const char *getMateReferenceSequence();
    virtual genomeIndex_t getMatePosition();
    virtual int getISize();     // insertion size -> 0 for single, delta for paired
    const char *getSequence(std::string &sequence, genomeIndex_t genomeMatchPosition);
    virtual bool hasMate();
    virtual MatchedReadBase & getMate();
#endif

};


#endif /* _MATCHEDREADBASE_H_ */
