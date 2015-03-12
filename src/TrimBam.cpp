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
#include "PhoneHome.h"
#include "SamFilter.h"

void TrimBam::trimBamDescription()
{
    std::cerr << " trimBam - Trim the ends of reads in a SAM/BAM file changing read ends to 'N' and quality to '!' or softclipping the ends (resulting file will not be sorted)" << std::endl;
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
    std::cerr << "By default reverse strands are reversed and then the left & right are trimmed .\n";
    std::cerr << "This means that --left actually trims from the right of the read in the SAM/BAM for reverse reads.\n";
    std::cerr << "Optionally --ignoreStrand/-i can be specified to ignore the strand information and treat forward/reverse the same.\n";
    std::cerr << "trimBam will modify the sequences to 'N', and the quality string to '!', unless --clip/-c is specified.\n";
    std::cerr << "--clip/-c indicates to soft clip instead of modifying the sequence or quality\n";
    std::cerr << "\tWhen clipping:\n";
    std::cerr << "\t  * if the entire read would be soft clipped, no clipping is done, and instead the read is marked as unmapped\n";
    std::cerr << "\t  * mate information is not updated (start positions/mapping may change after soft clipping)\n";
    std::cerr << "\t        * run samtools fixmate to fix mate information (will first need to sort by read name)\n";
    std::cerr << "\t  * output is not sorted (start positions/mapping may change after soft clipping)\n";
    std::cerr << "\t        * run samtools sort to resort by coordinate (after fixmate)\n";
    std::cerr << "\t  * soft clips already in the read are maintained or added to\n";
    std::cerr << "\t        * if 3 bases were clipped and 2 are specified to be clipped, no change is made to that end\n";
    std::cerr << "\t        * if 3 bases were clipped and 5 are specified to be clipped, 2 additional bases are clipped from that end\n";
}

// main function
int TrimBam::execute(int argc, char ** argv)
{
  SamFile samIn;
  SamFile samOut;
  int numTrimBaseL = 0;
  int numTrimBaseR = 0;
  bool noeof = false;
  bool ignoreStrand = false;
  bool clip = false;
  bool noPhoneHome = false;
  std::string inName = "";
  std::string outName = "";

  if ( argc < 5 ) {
    usage();
    std::cerr << "ERROR: Incorrect number of parameters specified\n";
    return(-1);
  }
  inName = argv[2];
  outName = argv[3];

  static struct option getopt_long_options[] = {
      // Input options
      { "left", required_argument, NULL, 'L'},
      { "right", required_argument, NULL, 'R'},
      { "ignoreStrand", no_argument, NULL, 'i'},
      { "clip", no_argument, NULL, 'c'},
      { "noeof", no_argument, NULL, 'n'},
      { "noPhoneHome", no_argument, NULL, 'p'},
      { "nophonehome", no_argument, NULL, 'P'},
      { "phoneHomeThinning", required_argument, NULL, 't'},
      { "phonehomethinning", required_argument, NULL, 'T'},
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
                            "L:R:icn", getopt_long_options, &n_option_index) )
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
          case 'i':
              ignoreStrand = true;
              break;
          case 'c':
              clip = true;
              break;
          case 'n':
              noeof = true;
              break;
          case 'p':
          case 'P':
              noPhoneHome = true;
              break;
          case 't':
          case 'T':
              PhoneHome::allThinning = atoi(optarg);
              break;
          default:
              fprintf(stderr,"ERROR: Unrecognized option %s\n",
                      getopt_long_options[n_option_index].name);
              return(-1);
      }
  }

  if(!noPhoneHome)
  {
      PhoneHome::checkVersion(getProgramName(), VERSION);
  }
  
  if(noeof)
  {
      // Set that the eof block is not required.
      BgzfFileType::setRequireEofBlock(false);
  }

  if ( ! samIn.OpenForRead(inName.c_str()) ) {
      fprintf(stderr, "***Problem opening %s\n",inName.c_str());
    return(-1);
  }

  if(!samOut.OpenForWrite(outName.c_str())) {
    fprintf(stderr, "%s\n", samOut.GetStatusMessage());
    return(samOut.GetStatus());
  }
  
  fprintf(stderr,"Arguments in effect: \n");
  fprintf(stderr,"\tInput file : %s\n",inName.c_str());
  fprintf(stderr,"\tOutput file : %s\n",outName.c_str());
  std::string trimType = "trim";
  if(clip) { trimType = "clip"; }
  if(numTrimBaseL == numTrimBaseR)
  {
      fprintf(stderr,"\t#Bases to %s from each side : %d\n", 
              trimType.c_str(), numTrimBaseL);
  }
  else
  {
      fprintf(stderr,"\t#Bases to %s from the left of forward strands : %d\n",
              trimType.c_str(), numTrimBaseL);
      fprintf(stderr,"\t#Bases to %s from the right of forward strands: %d\n",
              trimType.c_str(), numTrimBaseR);
      if(!ignoreStrand)
      {
          // By default, reverse strands are treated the opposite.
          fprintf(stderr,"\t#Bases to %s from the left of reverse strands : %d\n",
                  trimType.c_str(), numTrimBaseR);
          fprintf(stderr,"\t#Bases to %s from the right of reverse strands : %d\n",
                  trimType.c_str(), numTrimBaseL);
      }
      else
      {
          // ignore strand, treating forward & reverse strands the same
          fprintf(stderr,"\t#Bases to %s from the left of reverse strands : %d\n",
                  trimType.c_str(), numTrimBaseL);
          fprintf(stderr,"\t#Bases to %s from the right of reverse strands : %d\n",
                  trimType.c_str(), numTrimBaseR);
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
     // set based on ignoreStrand flag and strand info.
     int trimLeft = numTrimBaseL;
     int trimRight = numTrimBaseR;
     if(!ignoreStrand)
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
         return(-1);
       }

       if(clip)
       {
           SamFilter::softClip(samRecord, trimLeft, trimRight);
       }
       else
       {
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
     }

     if(!samOut.WriteRecord(samHeader, samRecord)) {
         // Failed to write a record.
       fprintf(stderr, "Failure in writing record %s\n", samOut.GetStatusMessage());
       return(-1);
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
