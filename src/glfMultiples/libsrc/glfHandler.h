////////////////////////////////////////////////////////////////////// 
// libsrc/glfHandler.h 
// (c) 2000-2010 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday June 16, 2010
// 
 
#ifndef __GLF_HANDLER_H__
#define __GLF_HANDLER_H__

#include "InputFile.h"
#include "StringBasics.h"

#pragma pack(push)
#pragma pack(1)

struct glfIndel
   {
   // Likelihood for the 1/1, 2/2 and 1/2
   unsigned char lk[3];

   // Allele lengths
   short length[2];

   unsigned char padding[3];
   };

struct glfEntry
   {
   /**  "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
   unsigned char refBase:4, recordType:4;

   /** offset of this record from the previous one, in bases */
   unsigned int offset;

   /** log10 minimum likelihood * 10 and the number of mapped reads */
   unsigned depth:24, minLLK:8;

   /** root mean squared maximum mapping quality for overlapping reads */
   unsigned char mapQuality;

   union {
   /** log10 likelihood ratio * 10 for genotypes AA, AC, AG, AT, CC, CG, CT, GG, GT, TT */
   unsigned char lk[10];
   glfIndel indel;
   };

   glfEntry & operator = (glfEntry & rhs);
   };

#pragma pack(pop)

class glfHandler
   {
   public:
      // Global information about the current GLF file
      IFILE    handle;
      String   header;

      // Information about the current section
      String   label;
      int      sections;
      int      currentSection;
      int      maxPosition;

      // Information on whether the end of the current section has been reached
      bool   endOfSection;

      // Currently active GLF record
      glfEntry data;
      int      position;
      double   likelihoods[10];
      String   indelSequence[2];

      // Error message in case previous command fails
      const char * errorMsg;

      glfHandler();
      ~glfHandler();

      bool Open(const String & filename);
      bool Create(const String & filename);
      bool isOpen();
      void Close();
      void Rewind();

      bool NextSection();
      bool NextEntry();
      bool NextBaseEntry();

      void BeginSection(const String & sectionLabel, int sectionLength);
      void EndSection();

      void WriteEntry(int outputPosition);

      char     GetReference(int position, char defaultBase);
      int      GetDepth(int position);
      const double * GetLikelihoods(int position);
      const unsigned char *   GetLogLikelihoods(int position);
      int      GetMapQuality(int position);

      static const double * GetDefaultLikelihoods()
         { return nullLikelihoods; }
      static const unsigned char * GetDefaultLogLikelihoods()
         { return nullLogLikelihoods; }

      static int GenotypeIndex(int base1, int base2)
         {
         return base1 < base2 ? (base1 - 1) * (10 - base1) / 2 + (base2 - base1) :
                                (base2 - 1) * (10 - base2) / 2 + (base1 - base2);
         }

   private:
      static char           translateBase[16];
      static char           backTranslateBase[5];
      static double         nullLikelihoods[10];
      static unsigned char  nullLogLikelihoods[10];

      bool ReadHeader();
      void WriteHeader(const String & headerText = "");
   };

#endif

 
