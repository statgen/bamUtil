////////////////////////////////////////////////////////////////////// 
// pdf/PDF.h 
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
 
#ifndef __PDF_H__
#define __PDF_H__

#include "PDFinfo.h"
#include "PDFpage.h"
#include "PDFfont.h"

#include <stdio.h>

class PDF
   {
   public:
      PDFPage page;
      PDFFont font;
      PDFInfo info;
      FILE *  file;

      PDF();
      ~PDF();

      void OpenFile(const char * name);
      void CloseFile();
      bool isOpen() { return file != NULL; }

      void WriteBoolean(bool boolean);
      void WriteInteger(int integer);
      void WriteDouble(double value);
      void WriteName(const char * name);
      void WriteString(const char * string);
      void WriteComment(const char * string);
      void WriteReference(int object);
      void WriteDate(int year, int month, int day);
      void WriteArray(const IntArray & array);
      void WriteReferenceArray(const IntArray & array);

      void WriteBoolean(const char * name, bool boolean);
      void WriteInteger(const char * name, int integer);
      void WriteDouble(const char * name, double value);
      void WriteName(const char * name, const char * name2);
      void WriteString(const char * name, const char * string);
      void WriteReference(const char * name, int object);
      void WriteDate(const char * name, int year, int month, int day);
      void WriteArray(const char * name, const IntArray & array);
      void WriteReferenceArray(const char * name, const IntArray & array);

      int  GetObject();
      void OpenObject(int object);
      void CloseObject();

      void OpenArray();
      void CloseArray();

      void OpenDictionary();
      void CloseDictionary();

      void WriteInteger(int object, int integer);
      void WriteString(int object, const char * string);

      int  OpenStream();
      void CloseStream();
      int  StreamLength();
      void AppendToStream(const char * string, ...);

      void LineBreak();

   private:
      IntArray objects;
      int length_index, stream_start;

      bool BalancedParenthesis(const char * string);
   };

#endif


 
 
 
 
