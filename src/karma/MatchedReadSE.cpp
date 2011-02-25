#include "MatchedReadSE.h"

void MatchedReadSE::updateMatch(MatchedReadSE& betterMatch)
{
    ((MatchedReadBase*)this)->updateMatch((MatchedReadBase&) betterMatch);
    return;
}

