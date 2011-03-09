#ifndef _MATCHEDREADSE_H_
#define _MATCHEDREADSE_H_

#include "MatchedReadBase.h"

class MatchedReadSE : public MatchedReadBase
{
 public:
    // overloaded methods - see base class for descriptions:
    // char *getSequence(String &sequence, genomeIndex_t genomeMatchPosition);
    void updateMatch(MatchedReadSE& betterMatch);

};


#endif /* _MATCHEDREADSE_H_ */
