////////////////////////////////////////////////////////////////////// 
// thunder/ShotgunManners.cpp 
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
 
#include "Manners.h"

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

static const char * machCrashExplanation = NULL;

void SetCrashExplanation(const char * specialMessage)
   {
   machCrashExplanation = specialMessage;
   }

void UnsetCrashExplanation()
   {
   machCrashExplanation = NULL;
   }

void SetupCrashHandlers()
   {
   signal(SIGINT, (signal_handler) UserBreak);
   signal(SIGSEGV, (signal_handler) OutOfMemory);
   signal(SIGABRT, (signal_handler) OutOfMemory);
   }

void OutOfMemory(int)
   {
   printf("\n\nMACH 1.0 FOR SHOTGUN SEQUENCES HAS CRASHED\n\n"
          "The operating system has decided to terminate this run,\n"
          "of the Markov Chain Haplotyper (MACH). Either MACH 1.0\n"
          "requested too much memory or you have encountered a bug.\n\n"
          "If you don't think this is a memory issue, you can help\n"
          "improve this program by reporting bugs via a short e-mail\n"
          "to goncalo@umich.edu. These e-mails are most helpful if you\n"
          "include a description of the problem and example of how it can\n"
          "be reproduced.\n\n");

   if (machCrashExplanation != NULL)
      printf("MACH 1.0 crashed while %s\n", machCrashExplanation);

   exit(EXIT_FAILURE);
   }

void UserBreak(int)
   {
   printf("\n\nMACH 1.0 STOPPED BY USER\n\n");

   exit(EXIT_FAILURE);
   }

 
