#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include "StringArray.h"

class Sequence
{
 public:
  static char complement(char);
  static int nt2idx(char);
  static int nt2int(char);
  static String reverse(String &S);
  static String revComp(String &S);
};

#endif
