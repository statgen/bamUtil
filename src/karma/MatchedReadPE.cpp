#include <iostream>
#include <string>
#include "MatchedReadPE.h"

//
//
void MatchedReadPE::constructorClear()
{
    MatchedReadBase::constructorClear();

    pairCumulativePosteriorProbabilities = 0.0;
    pairQuality = UNSET_QUALITY;

}

double MatchedReadPE::getQualityScore()
{
    return MatchedReadBase::getQualityScore(pairQuality, pairCumulativePosteriorProbabilities);
}


void MatchedReadPE::printOptionalTags(std::ostream &file, bool isPairAligned, const std::string &sampleGroupID, const std::string &alignmentPathTag)
{
    MatchedReadBase::printOptionalTags(file, isPairAligned, sampleGroupID, alignmentPathTag);
    if (isPairAligned && pairQuality) file << "\tPQ:i:" <<  pairQuality;
}

///
/// @param betterMatch a MatchedReadPE value used to update this class
/// Update our bestMatch object with the newer information.
/// Called while we are iterating over all possible matches,
/// and we discover what we think is a better one.
///
/// We could in thoery copy over nearly all the information,
/// but most is not used - just the pairQuality, genome match
/// position and soon a few others if we do posterior quality
/// checks.
///
void MatchedReadPE::updateMatch(MatchedReadPE & betterMatch)
{
    MatchedReadBase::updateMatch(betterMatch);
    pairQuality = betterMatch.pairQuality;
}

