/*
 *  Copyright (C) 2010-2012  Regents of the University of Michigan,
 *                           Hyun Min Kang
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

// Modified 4/3/12 by Mary Kate Trost to put into bamUtil.

#include <stdio.h>
#include <string.h>
#include "SamFile.h"
#include "TrimBam.h"

void TrimBam::trimBamDescription()
{
    std::cerr << " trimBam - Trim the ends of reads in a SAM/BAM file changing read ends to 'N' and quality to '!'" << std::endl;
}


void TrimBam::description()
{
    trimBamDescription();
}


void TrimBam::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam trimBam [inFile] [outFile] [num-bases-to-trim-on-each-side]\n";
    std::cerr << "trimBam will modify the sequences to 'N', and the quality string to '!'\n";
}

// main function
int TrimBam::execute(int argc, char ** argv)
{
  SamFile samIn;
  SamFile samOut;
  int numTrimBase;

  if ( argc != 5 ) {
    usage();
    abort();
  }

  if ( ! samIn.OpenForRead(argv[2]) ) {
      fprintf(stderr, "***Problem opening %s\n",argv[1]);
    abort();
  }

  if(!samOut.OpenForWrite(argv[3])) {
    fprintf(stderr, "%s\n", samOut.GetStatusMessage());
    return(samOut.GetStatus());
  }
  
  numTrimBase = atoi(argv[4]);

  fprintf(stderr,"Arguments in effect: \n");
  fprintf(stderr,"\tInput file : %s\n",argv[2]);
  fprintf(stderr,"\tOutput file : %s\n",argv[3]);
  fprintf(stderr,"\t#TrimBases : %d\n",numTrimBase);
  
   // Read the sam header.
   SamFileHeader samHeader;
   if(!samIn.ReadHeader(samHeader))
   {
      fprintf(stderr, "%s\n", samIn.GetStatusMessage());
      return(samIn.GetStatus());
   }

   // Write the sam header.
   if(!samOut.WriteHeader(samHeader))
   {
      fprintf(stderr, "%s\n", samOut.GetStatusMessage());
      return(samOut.GetStatus());     
   }

   SamRecord samRecord;
   char seq[65536];
   char qual[65536];
   int i, len;

   // Keep reading records until ReadRecord returns false.
   while(samIn.ReadRecord(samHeader, samRecord)) {
     // Successfully read a record from the file, so write it.
     strcpy(seq,samRecord.getSequence());
     strcpy(qual,samRecord.getQuality());
     len = strlen(seq);
     // Do not trim if sequence is '*'
     if ( strcmp(seq, "*") != 0 ) {
       bool qualValue = true;
       if(strcmp(qual, "*") == 0)
       {
           qualValue = false;
       }
       int qualLen = strlen(qual);
       if ( (qualLen != len) && qualValue ) {
         fprintf(stderr,"ERROR: Sequence and Quality have different length\n");
         abort();
       }
       if ( len < numTrimBase ) {
         for(i=0; i < len; ++i) {
           seq[i] = 'N';
           if ( qualValue ) {
             qual[i] = '!';
           }
         }
       }
       else {
         for(i=0; i < numTrimBase; ++i) {
           seq[i] = 'N';
           seq[len-i-1] = 'N';
           if ( qualValue ) {
             qual[i] = '!';
             qual[len-i-1] = '!';
           }
         }
       }
       samRecord.setSequence(seq);
       samRecord.setQuality(qual);
     }

     if(!samOut.WriteRecord(samHeader, samRecord)) {
         // Failed to write a record.
       fprintf(stderr, "Failure in writing record %s\n", samOut.GetStatusMessage());
       abort();
     }
   }
   
   if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
   {
      // Failed to read a record.
      fprintf(stderr, "%s\n", samIn.GetStatusMessage());
   }   
   
   std::cerr << std::endl << "Number of records read = " << 
     samIn.GetCurrentRecordCount() << std::endl;
   std::cerr << "Number of records written = " << 
     samOut.GetCurrentRecordCount() << std::endl;

   if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
   {
     // Failed reading a record.
     return(samIn.GetStatus());
   }

   // Since the reads were successful, return the status based
   samIn.Close();
   samOut.Close();
   return 0;
}
