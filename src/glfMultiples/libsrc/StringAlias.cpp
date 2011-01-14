////////////////////////////////////////////////////////////////////// 
// libsrc/StringAlias.cpp 
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
 
#include "StringAlias.h"
#include "InputFile.h"

void StringAlias::SetAlias(String & string, String & alias)
   {
   int index = lookup.Integer(string);

   if (index < 0)
      {
      aliases.Push(alias);
      lookup.SetInteger(string, aliases.Length() - 1);
      }
   else
      aliases[index] = alias;
   }

const String & StringAlias::GetAlias(const String & string) const
   {
   int index = lookup.Integer(string);

   if (index < 0)
      return string;
   else
      return aliases[index];
   }

bool StringAlias::ReadFromFile(const char * filename)
   {
   IFILE input = ifopen(filename, "rt");

   if (input == NULL)
      return false;

   ReadFromFile(input);

   ifclose(input);

   return true;
   }

bool StringAlias::ReadFromFile(IFILE & input)
   {
   StringArray lines, tokens;
   lines.Read(input);

   for (int j = 0; j < lines.Length(); j++)
      {
      tokens.ReplaceTokens(lines[j]);

      if (tokens.Length() != 2) continue;

      SetAlias(tokens[0], tokens[1]);
      }

   return true;
   }

 
