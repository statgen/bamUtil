////////////////////////////////////////////////////////////////////// 
// pdf/PDFinfo.cpp 
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
 
 
#include "PDFinfo.h"
#include "PDF.h"

PDFInfo::PDFInfo()
   {
   time_t creation = time(NULL);
   tm * time_struct = gmtime(&creation);

   creationYear = time_struct->tm_year + 1900;
   creationMonth = time_struct->tm_mon + 1;
   creationDay = time_struct->tm_mday;

   creator = "A C++ PDF Library (c) 2001 Goncalo Abecasis";
   }

void PDFInfo::Write(PDF & pdf)
   {
   object = pdf.GetObject();

   pdf.OpenObject(object);
   pdf.OpenDictionary();

   if (!title.IsEmpty()) pdf.WriteString("Title", title);
   if (!author.IsEmpty()) pdf.WriteString("Author", author);
   if (!subject.IsEmpty()) pdf.WriteString("Subject", subject);
   if (!keywords.IsEmpty()) pdf.WriteString("Keywords", keywords);
   if (!producer.IsEmpty()) pdf.WriteString("Producer", producer);

   pdf.WriteString("Creator", creator);
   pdf.WriteDate("CreationDate", creationYear, creationMonth, creationDay);

   pdf.CloseDictionary();
   pdf.CloseObject();
   } 
 
 
 
