#include "QSamFlag.h"


QSamFlag::QSamFlag()
{
  isRead1 = isRead2 = isPaired = isUnPaired = isProperPaired = isUnMapped = isReverse = isDup = isQCFail = false;
}
void QSamFlag::GetFlagFields(int flag)
{
  if(flag&BAM_FPAIRED) 
    {
      isPaired = true;
      if(flag & BAM_FREAD1) isRead1 = true;
      if(flag & BAM_FREAD2) isRead2 = true;
      if(flag & BAM_FPROPER_PAIR) isProperPaired = true;
    }

  if(flag & BAM_FUNMAP) isUnMapped = true;
  if(flag & BAM_FREVERSE) isReverse = true;
  if(flag & BAM_FDUP) isDup = true;
  if(flag & BAM_FQCFAIL) isQCFail = true;
}
