////////////////////////////////////////////////////////////////////// 
// pdf/PDFinfo.h 
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
 
#ifndef __PDFINFO_H__
#define __PDFINFO_H__

#include "StringBasics.h"

#include <time.h>

class PDF;

class PDFInfo
   {
   public:
      String title;
      String author;
      String subject;
      String keywords;
      String producer;
      int    object;

      PDFInfo();

      void Write(PDF & pdf);

   private:
      String creator;
      int    creationYear, creationMonth, creationDay;
   };

#endif

 
 
 
 
