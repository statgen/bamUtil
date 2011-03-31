////////////////////////////////////////////////////////////////////// 
// mach1/ErrorRate.h 
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
 
#ifndef __ERRORRATE_H__
#define __ERRORRATE_H__

class Errors
   {
   public:
      int    mismatches;
      int    matches;
      int    uncertain_pairs;

      float  rate;

      Errors();

      float  Update();
      void   Reset();

      Errors & operator += (const Errors & rhs);
   };

#endif

 
