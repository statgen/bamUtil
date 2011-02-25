#ifndef _MATCHEDREADPE_H_
#define _MATCHEDREADPE_H_

#include "MatchedReadBase.h"

//
// these are actually match candidates, and serve three
// purposes:
// 1) running bestMatch
// 2) working currentMatch
// 3) part of a list of possible matches, which we
//    evaluate the Q score for in a lazy fashion.
//    In this instance, we need all args to pass to wordIndex::getQvalue
//
class MatchedReadPE : public MatchedReadBase
{
public:
    int pairQuality;
    inline bool pairQualityIsValid() const
    {
        return pairQuality >=0;
    }
    double pairCumulativePosteriorProbabilities;

    MatchedReadPE()
    {
        constructorClear();
    }
    void constructorClear();
    void printOptionalTags(
        std::ostream &,
        bool isProperAligned,
        const std::string &sampleGroupID,
        const std::string &alignmentPathTag);
    void updateMatch(MatchedReadPE &);
    double getQualityScore();
};



#endif /* _MATCHEDREADPE_H_ */
