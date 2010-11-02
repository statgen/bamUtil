#include "Sequence.h"

char Sequence::complement(char c)
{
 switch(toupper(c)){
        case 'A': return('T');
        case 'C': return('G');
        case 'G': return('C');
        case 'T': return('A');
        default:  return(toupper(c));
        }
}

String Sequence::reverse(String & s)
{
  String tmp = s;
  int L = s.Length();
  for(int i=0; i<L; i++)
    tmp[i] = s[L-i-1];
  return(tmp);
}

String Sequence::revComp(String & s)
{
  String tmp = s;
  int L = s.Length();
  for(int i=0; i<L; i++)
    tmp[i] = complement(s[L-i-1]);
  return(tmp);
}

int Sequence::nt2idx(char c) {
  switch(toupper(c)){
  case 'A': return(0);
  case 'C': return(1);
  case 'G': return(2);
  case 'T': return(3);
  case 'N': return(4);
  default:  return(-1);
  }
}
                                                                

int Sequence::nt2int(char c) {
  switch(toupper(c)){
  case 'A': return(1);
  case 'C': return(2);
  case 'G': return(3);
  case 'T': return(4);
  case 'N': return(5);
  default:  return(-1);
  }
}
                                                                

