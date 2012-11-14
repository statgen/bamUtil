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
#include <getopt.h>
#include "SamFile.h"
#include "SamFlag.h"
#include "BgzfFileType.h"
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
    std::cerr << "Alternately, the number of bases from each side can be specified (either or both -L/-R (--left/--right) can be specified):\n";
    std::cerr << "\t./bam trimBam [inFile] [outFile] -L [num-bases-to-trim-from-left] -R [num-bases-to-trim-from-right]\n";
    std::cerr << "By default Left/Right is as the reads are in the SAM/BAM file.\n";
    std::cerr << "Optionally --reverse/-r can be specified to reverse the left/right for reverse reads\n";
    std::cerr << "trimBam will modify the sequences to 'N', and the quality string to '!'\n";
}

// main function
int TrimBam::execute(int argc, char ** argv)
{
  SamFile samIn;
  SamFile samOut;
  int numTrimBaseL = 0;
  int numTrimBaseR = 0;
  bool noeof = false;
  bool reverse = false;
  std::string inName = "";
  std::string outName = "";

  if ( argc < 5 ) {
    usage();
    abort();
  }
  inName = argv[2];
  outName = argv[3];

  static struct option getopt_long_options[] = {
      // Input options
      { "left", required_argument, NULL, 'L'},
      { "right", required_argument, NULL, 'R'},
      { "reverse", no_argument, NULL, 'r'},
      { "noeof", no_argument, NULL, 'n'},
      { NULL, 0, NULL, 0 },
  };
  
  int argIndex = 4;
  if(argv[argIndex][0] != '-')
  {
      // This is the number of bases to trim off both sides
      // so convert to a number.
      numTrimBaseL = atoi(argv[argIndex]);
      numTrimBaseR = numTrimBaseL;
      ++argIndex;
  }

  int c = 0;
 int n_option_index = 0;
  // Process any additional parameters
  while ( ( c = getopt_long(argc, argv,
                            "L:R:rn", getopt_long_options, &n_option_index) )
          != -1 )
  {
      switch(c) 
      {
          case 'L':
              numTrimBaseL = atoi(optarg);
              break;
          case 'R':
              numTrimBaseR = atoi(optarg);
              break;
          case 'r':
              reverse = true;
              break;
          case 'n':
              noeof = true;
              break;
          default:
              fprintf(stderr,"Unrecognized option %s",
                      getopt_long_options[n_option_index].name);
              abort();
      }
  }

  if(noeof)
  {
      // Set that the eof block is not required.
      BgzfFileType::setRequireEofBlock(false);
  }

  if ( ! samIn.OpenForRead(inName.c_str()) ) {
      fprintf(stderr, "***Problem opening %s\n",inName.c_str());
    abort();
  }

  if(!samOut.OpenForWrite(outName.c_str())) {
    fprintf(stderr, "%s\n", samOut.GetStatusMessage());
    return(samOut.GetStatus());
  }
  
  fprintf(stderr,"Arguments in effect: \n");
  fprintf(stderr,"\tInput file : %s\n",inName.c_str());
  fprintf(stderr,"\tOutput file : %s\n",outName.c_str());
  if(numTrimBaseL == numTrimBaseR)
  {
      fprintf(stderr,"\t#Bases to trim from each side : %d\n", numTrimBaseL);
  }
  else
  {
      fprintf(stderr,"\t#Bases to trim from the left of forward strands : %d\n",
              numTrimBaseL);
      fprintf(stderr,"\t#Bases to trim from the right of forward strands: %d\n",
              numTrimBaseR);
      if(reverse)
      {
          // When using the reverse option, reverse strands are treated the opposite.
          fprintf(stderr,"\t#Bases to trim from the left of reverse strands : %d\n",
                  numTrimBaseR);
          fprintf(stderr,"\t#Bases to trim from the right of reverse strands : %d\n",
                  numTrimBaseL);
      }
      else
      {
          fprintf(stderr,"\t#Bases to trim from the left of reverse strands : %d\n",
                  numTrimBaseL);
          fprintf(stderr,"\t#Bases to trim from the right of reverse strands : %d\n",
                  numTrimBaseR);
      }
  }
 
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

     // Number of bases to trim from the left/right,
     // set based on reverse flag and strand info.
     int trimLeft = numTrimBaseL;
     int trimRight = numTrimBaseR;
     if(reverse)
     {
         if(SamFlag::isReverse(samRecord.getFlag()))
         {
             // We are reversing the reverse reads,
             // so swap the left & right trim counts.
             trimRight = numTrimBaseL;
             trimLeft = numTrimBaseR;
         }
     }

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
       if ( len < (trimLeft + trimRight) ) {
         // Read Length is less than the total number of bases to trim,
         // so trim the entire read.
         for(i=0; i < len; ++i) {
           seq[i] = 'N';
           if ( qualValue ) {
             qual[i] = '!';
           }
         }
       }
       else
       {
           // Read Length is larger than the total number of bases to trim,
           // so trim from the left, then from the right.
           for(i=0; i < trimLeft; ++i)
           {
               // Trim the bases from the left.
               seq[i] = 'N';
               if ( qualValue )
               {
                   qual[i] = '!';
               }
           }
           for(i = 0; i < trimRight; i++)
           {
               seq[len-i-1] = 'N';
               if(qualValue)
               {
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
