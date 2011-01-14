////////////////////////////////////////////////////////////////////// 
// libsrc/BaseQualityHelper.cpp 
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
 
#include "BaseQualityHelper.h"

#include <math.h>

baseQualityConvertor bQualityConvertor;

baseQualityConvertor::baseQualityConvertor()
   {
   // Create a quick lookup table to speed up conversion of
   // base quality values stored as log10 (error rates) into
   // fractional error rates
   for (int i = 0; i <= 255; i++)
      doubleLookup[i] = pow(0.1, i * 0.1);
   // doubleLookup[255] = 0.0;
   }

double baseQualityConvertor::toDouble(unsigned char bq)
   {
   return doubleLookup[bq];
   }

 
