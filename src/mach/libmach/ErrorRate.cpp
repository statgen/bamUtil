////////////////////////////////////////////////////////////////////// 
// mach1/ErrorRate.cpp 
// (c) 2000-2008 Goncalo Abecasis
// 
// This file is distributed as part of the MaCH source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile MaCH.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Saturday April 12, 2008
// 
 
#include "ErrorRate.h"

#include <math.h>

Errors::Errors()
   {
   rate = 0.0;

   Reset();
   }

float Errors::Update()
   {
   if (matches + mismatches > 0)
      {
      float previous = 0.0;
      rate = mismatches / (double) (matches + mismatches);

      if (uncertain_pairs)
         while ((rate > 1e-10) && (fabs(rate - previous) > rate * 1e-4))
            {
            double ratio = rate * rate / (rate * rate + (1.0 - rate) * (1.0 - rate));

            previous = rate;
            rate = (mismatches + ratio * uncertain_pairs * 2.0) / (matches + mismatches + uncertain_pairs * 2.0);
            }
      }
   else if (uncertain_pairs)
      rate = 0.0;

   return rate;
   }

void Errors::Reset()
   {
   matches = mismatches = uncertain_pairs = 0;
   }

Errors & Errors::operator += (const Errors & rhs)
   {
   matches += rhs.matches;
   mismatches += rhs.mismatches;
   uncertain_pairs += rhs.uncertain_pairs;

   return *this;
   }

 
