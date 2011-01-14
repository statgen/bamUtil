////////////////////////////////////////////////////////////////////// 
// libsrc/StringAlias.h 
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
 
#ifndef __STRINGALIAS_H__
#define __STRINGALIAS_H__

#include "StringArray.h"
#include "StringHash.h"

class StringAlias
   {
   public:
      StringAlias()               {}
      virtual ~StringAlias()      {}

      void SetAlias(String & string, String & alias);

      const String & GetAlias(const String & string) const;

      bool  ReadFromFile(const char * filename);
      bool  ReadFromFile(IFILE & input);

   private:
      StringIntHash  lookup;
      StringArray    aliases;
   };

#endif

 
