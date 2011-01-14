////////////////////////////////////////////////////////////////////// 
// pdf/PDF.cpp 
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
 
#include "PDF.h"
#include "Error.h"

#include <stdarg.h>

PDF::PDF() : page(*this), font(*this)
   {
   file = NULL;
   objects.Push(0);
   }

PDF::~PDF()
   {
   CloseFile();
   }

void PDF::OpenFile(const char * name)
   {
   file = fopen(name, "wb");

   if (file == NULL)
      error("Error opening file [%s]\n", name);

   char signature[] = {'%', '%', (char) ('G' + 128), (char) ('R' + 128), (char) ('A' + 128), '\n', '\n', 0};

   fprintf(file, "%%PDF-1.4\n");
   fprintf(file, signature);
   }

void PDF::CloseFile()
   {
   if (file == NULL)
      return;

   page.ClosePage();

   // If the PDF document is empty, add an empty
   // page to ensure well formed document.
   if (page.GetPageNumber() == 0)
      {
      page.OpenPage();
      page.ClosePage();
      }

   info.Write(*this);
   page.WritePageTree();
   font.WriteDictionary();

   int catalogue = GetObject();
   OpenObject(catalogue);
   OpenDictionary();

   WriteName("Type", "Catalog");
   WriteReference("Pages", page.tree_index);

   CloseDictionary();
   CloseObject();

   int trailer_offset = ftell(file);

   fprintf(file, "xref\n");
   fprintf(file, "0 %d\n", objects.Length() );
   fprintf(file, "0000000000 65535 f \n");

   for (int i = 1; i < objects.Length(); i++)
      fprintf(file, "%010d 00000 n \n", objects[i]);

   fprintf(file, "trailer\n");

   fprintf(file, "<<\n");
   fprintf(file, "/Size %d\n", objects.Length());
   fprintf(file, "/Root %d 0 R\n", catalogue);
   fprintf(file, "/Info %d 0 R\n", info.object);
   fprintf(file, ">>\n");

   fprintf(file, "startxref\n");
   fprintf(file, "%d\n", trailer_offset);
   fprintf(file, "%%%%EOF");

   fclose(file);
   file = NULL;
   }

bool PDF::BalancedParenthesis(const char * string)
   {
   int balance = 0;

   while (true)
      {
      switch (*string)
         {
         case '(' :
            balance++;
            break;
         case ')' :
            balance--;
            if (balance < 0)
               return false;
            break;
         case 0 :
            return balance == 0;
         default:
            break;
         }
      string++;
      }
   }

void PDF::WriteBoolean(bool boolean)
   {
   fprintf(file, " %s ", boolean ? "true" : "false");
   }

void PDF::WriteInteger(int integer)
   {
   fprintf(file, " %d ", integer);
   }

void PDF::WriteDouble(double value)
   {
   fprintf(file, " %#f ", value);
   }

void PDF::WriteName(const char * name)
   {
   fputc('/', file);
   while (*name)
      {
      if (*name >= 33 && *name <= 126 && *name != '#')
         fputc(*name, file);
      else
         fprintf(file, "#%02X", *name);
      name++;
      }
   fputc(' ', file);
   }

void PDF::WriteString(const char * string)
   {
   fputc('(', file);

   bool balanced = false, check = false; // Initialization avoids compiler warning

   while (*string)
      {
      switch (*string)
         {
         case '(' :
         case ')' :
            if (!check)
               {
               check = true;
               balanced = BalancedParenthesis(string);
               }
            if (balanced)
               {
               fputc(*string, file);
               break;
               }
         case '\\' :
            fputc('\\', file);
         default :
            fputc(*string, file);
         }
      string++;
      }

   fputc(')', file);
   fputc(' ', file);
   }

void PDF::WriteComment(const char * comment)
   {
   fputc('%', file);

   while (*comment)
      {
      int ch = *comment++;

      if (ch == '\r') continue;
      fputc(ch, file);
      if (ch == '\n') fputc('%', file);
      }
   }

void PDF::OpenArray()
   {
   fprintf(file, " [ ");
   }

void PDF::CloseArray()
   {
   fprintf(file, " ]\n");
   }

void PDF::WriteArray(const IntArray & array)
   {
   fprintf(file, " [ ");

   for (int i = 0; i < array.Length(); i++)
      WriteInteger(array[i]);

   fprintf(file, " ] ");
   }

void PDF::WriteReferenceArray(const IntArray & array)
   {
   fprintf(file, " [ ");

   for (int i = 0; i < array.Length(); i++)
      WriteReference(array[i]);

   fprintf(file, " ] ");
   }

void PDF::WriteReference(int object)
   {
   fprintf(file, " %d 0 R ", object);
   }

void PDF::WriteDate(int year, int month, int day)
   {
   fprintf(file, " (D:%04d%02d%02d) ", year, month, day);
   }

void PDF::WriteBoolean(const char * name, bool boolean)
   {
   WriteName(name);
   WriteBoolean(boolean);
   fputc('\n', file);
   }

void PDF::WriteInteger(const char * name, int integer)
   {
   WriteName(name);
   WriteInteger(integer);
   fputc('\n', file);
   }

void PDF::WriteDouble(const char * name, double value)
   {
   WriteName(name);
   WriteDouble(value);
   fputc('\n', file);
   }

void PDF::WriteName(const char * name, const char * name2)
   {
   WriteName(name);
   WriteName(name2);
   fputc('\n', file);
   }

void PDF::WriteString(const char * name, const char * string)
   {
   WriteName(name);
   WriteString(string);
   fputc('\n', file);
   }

void PDF::WriteReference(const char * name, int object)
   {
   WriteName(name);
   WriteReference(object);
   fputc('\n', file);
   }

void PDF::WriteDate(const char * name, int year, int month, int day)
   {
   WriteName(name);
   WriteDate(year, month, day);
   fputc('\n', file);
   }

void PDF::WriteArray(const char * name, const IntArray & array)
   {
   WriteName(name);
   WriteArray(array);
   fputc('\n', file);
   }

void PDF::WriteReferenceArray(const char * name, const IntArray & array)
   {
   WriteName(name);
   WriteReferenceArray(array);
   fputc('\n', file);
   }


int PDF::GetObject()
   {
   objects.Push(0);
   return objects.Length() - 1;
   };

void PDF::OpenObject(int object)
   {
   objects[object] = ftell(file);

   fprintf(file, "%d 0 obj\n", object);
   };

void PDF::CloseObject()
   {
   fprintf(file, "\nendobj\n\n");
   }

void PDF::WriteInteger(int object, int integer)
   {
   OpenObject(object);
   WriteInteger(integer);
   CloseObject();
   }

void PDF::WriteString(int object, const char * string)
   {
   OpenObject(object);
   WriteString(string);
   CloseObject();
   }

int PDF::OpenStream()
   {
   int object = GetObject();
   length_index = GetObject();

   OpenObject(object);
   fprintf(file, " << /Length %d 0 R >>\n", length_index);

   fprintf(file, "stream\n");
   stream_start = ftell(file);

   return object;
   }

void PDF::CloseStream()
   {
   int length = ftell(file) - stream_start;

   fprintf(file, "\nendstream");
   CloseObject();

   WriteInteger(length_index, length);
   }

int PDF::StreamLength()
   {
   return ftell(file) - stream_start;
   }

void PDF::AppendToStream(const char * string, ...)
   {
   va_list argptr;

   va_start(argptr, string);
   vfprintf(file, string, argptr);
   va_end(argptr);
   }

void PDF::OpenDictionary()
   {
   fprintf(file, "<<\n");
   }

void PDF::CloseDictionary()
   {
   fprintf(file, ">>");
   }

void PDF::LineBreak()
   {
   fprintf(file, "\n");
   }
 
 
 
 
