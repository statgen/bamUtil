#ifndef __SAMFLAG_H__
#define __SAMFLAG_H__

#include "FlagDef.h"

class QSamFlag
{
public:
  bool isRead1;
    bool isRead2;
      bool isPaired;
      bool isUnPaired;
      bool isProperPaired;
        bool isUnMapped;
          bool isReverse;
            bool isDup;
              bool isQCFail;
              
              public:
                QSamFlag();
                  void GetFlagFields(int);
                  void SetRead1(bool b){ isRead1 = b;};
                  void SetRead2(bool b){ isRead2 = b;};
                  void SetUnMapped(bool b){ isUnMapped = b;};
                  void SetPaired(bool b){ isPaired = b; };
                  void SetUnPaired(bool b){ isUnPaired = b; };
                  void SetReverse(bool b){ isReverse = b;};
                  void SetDuplicate(bool b){ isDup = b;};
                  void SetQCFail(bool b) { isQCFail = b;};
                  };
                  

#endif
