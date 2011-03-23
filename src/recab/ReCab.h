/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
////////////////////////////////////////////////////////////////////////
// Recalibrator
//
// by Christian Fuchsberger
//
// - Last modified on June 9th, 2010
//
////////////////////////////////////////////////////////////////////////

#ifndef __RE_CAP_H__
#define __RE_CAP_H__


#include <stdint.h>
#include <string>
// imports from samtools
#include "SamFile.h"
#include "Generic.h"
#include "GenomeSequence.h"
#include "MemoryMapArray.h"
#include "HashErrorModel.h"
#include "Prediction.h"
#include "BaseAsciiMap.h"


class ReCab {
public:
	// quality String
	typedef struct {
	    std::string oldq;
	    std::string newq;
	} quality_t;

    //quality fields
	std::string qField;


   //stats
   uint64_t basecounts;
   uint64_t mappedCount;
   uint64_t unMappedCount;
   uint64_t mappedCountQ;

   uint64_t BunMappedCount;
   uint64_t BMappedCount;

   uint64_t  zeroMapQualCount;

   GenomeSequence referenceGenome;
   mmapArrayBool_t dbSNP;
   HashErrorModel hasherrormodel;
   Prediction prediction;

   ReCab();
   ~ReCab();
   // conversion table
   static int nt2idx2[256];
   static void conversionTable();
   static int nt2idx(char c);
   static char complement(char c);
   inline bool processRead(SamRecord& record,int processtype,quality_t& quality_strings);
   static uint32_t addTokenizedStrings(const std::string& str, const std::string& delimiters, std::vector<std::string>& tokens);

private:
    BaseAsciiMap myBaseAsciiMap;
};

class CigarCol
{
 public:
  int *fieldCount;
  int size;
  String fieldTypes;
  String tmp;
  int nClips_begin, nClips_end;

  CigarCol() { fieldCount=NULL; nClips_begin=nClips_end=0;};
  ~CigarCol(){ CigarClear(); };

  void GetClipLength_begin(String cigar)
  {
    for(int i=0; i<cigar.Length(); i++)
    {
      if(cigar[i]!='S' && !isdigit(cigar[i])) return;
      if(cigar[i]=='S')
        {
          nClips_begin = cigar.SubStr(0, i).AsInteger();
          return;
        }
    }
 }
 void GetClipLength_end(String cigar)
{
	if(cigar[cigar.Length()-1] != 'S') return;

	for(int i=cigar.Length()-2; i>=0; i--)
	  if(!isdigit(cigar[i])) {
	  nClips_end = cigar.SubStr(i+1, cigar.Length()-i-2).AsInteger();
	  return; }
  }

  void CigarParse(String cigar) {
    size = 0;
    for  (int i = 0; i < cigar.Length(); i++)
      if (!isdigit(cigar[i])) size ++;

    fieldCount = new int[size];

    int k = 0; tmp = "";

    for  (int i = 0; i < cigar.Length(); i++)
      { if (isdigit(cigar[i])) tmp += cigar[i];
        else {
          fieldCount[k++] = tmp.AsInteger();
          tmp.Clear();
          fieldTypes += cigar[i];}
      }
  };

  void CigarClear(){
    if (fieldCount != NULL) {
      delete [] fieldCount;
      fieldCount = NULL;
      fieldTypes="";
      tmp="";}
  };
};



#endif
