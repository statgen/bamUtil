/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "StringArray.h"
#include "StringHash.h"
#include "Parameters.h"
#include "FastQFile.h"

int main(int argc, char ** argv)
{   
   ParameterList inputParameters;
   String filename;
   int minReadLength = 10;
   int printableErrors = 20;
   int maxErrors = -1;
   String testParam;
   BaseAsciiMap::SPACE_TYPE myBaseType = BaseAsciiMap::UNKNOWN;
   
   // Read the parameters from the command line.
   bool baseSpace = false;
   bool colorSpace = false;
   bool autoDetect = false;
   bool ignoreErrors = false;
   bool baseComposition = false;
   bool quiet = false;
   bool params = false;
   bool disableSeqIDCheck = false;

   BEGIN_LONG_PARAMETERS(longParameterList)
      LONG_STRINGPARAMETER("file", &filename)
      LONG_PARAMETER("baseComposition", &baseComposition)
      LONG_PARAMETER("disableSeqIDCheck", &disableSeqIDCheck)
      LONG_PARAMETER("quiet", &quiet)
      LONG_PARAMETER("params", &params)
      LONG_INTPARAMETER("minReadLen", &minReadLength)
      LONG_INTPARAMETER("maxErrors", &maxErrors)
      LONG_PARAMETER_GROUP("Space Type")
         EXCLUSIVE_PARAMETER("baseSpace", &baseSpace)
         EXCLUSIVE_PARAMETER("colorSpace", &colorSpace)
         EXCLUSIVE_PARAMETER("auto", &autoDetect)
      LONG_PARAMETER_GROUP("Errors")
         EXCLUSIVE_PARAMETER("ignoreErrors", &ignoreErrors)
         LONG_SMARTINTPARAMETER("printableErrors", &printableErrors)
   BEGIN_LEGACY_PARAMETERS()
      LONG_PARAMETER("printBaseComp", &baseComposition)       
      LONG_PARAMETER("disableAllMessages", &quiet)
      LONG_INTPARAMETER("quitAfterErrorNum", &maxErrors)
      LONG_PARAMETER_GROUP("Space Type")
         EXCLUSIVE_PARAMETER("baseSpace", &baseSpace)
         EXCLUSIVE_PARAMETER("colorSpace", &colorSpace)
         EXCLUSIVE_PARAMETER("autoDetect", &autoDetect)
      LONG_PARAMETER_GROUP("Errors")
         EXCLUSIVE_PARAMETER("ignoreAllErrors", &ignoreErrors)
         LONG_SMARTINTPARAMETER("maxReportedErrors", &printableErrors)
   END_LONG_PARAMETERS();
   
   inputParameters.Add(new LongParameters ("Input Parameters", longParameterList));

   inputParameters.Read(argc, argv);

   if(ignoreErrors)
   {
      // Ignore all errors, so set printableErrors to 0.
      printableErrors = 0;
   }

   // Set the base type based on the passed in parameters.
   if(baseSpace)
   {
      // Base Space
      myBaseType = BaseAsciiMap::BASE_SPACE;
   }
   else if(colorSpace)
   {
      myBaseType = BaseAsciiMap::COLOR_SPACE;
   }
   else
   {
      myBaseType = BaseAsciiMap::UNKNOWN;
      // Set autoDetect
      autoDetect = true;
   }

   // DO not print status if set to quiet.
   if((!quiet) && params)
   {
      inputParameters.Status();
   }

   if(filename == "")
   {
      if(quiet)
      {
         return(-1);
      }
      // No filename was specified so print a usage description.
      std::cout << "ERROR: No filename specified.  See below for usage help.";
      std::cout << std::endl << std::endl;

      std::cout << "  Required Parameters:" << std::endl;
      std::cout << "\t--file  :  FastQ filename with path to be prorcessed.\n";
      std::cout << std::endl;

      std::cout << "  Optional Parameters:" << std::endl;
      std::cout << "\t--minReadLen         : Minimum allowed read length (Defaults to 10).\n";
      std::cout << "\t--maxErrors          : Number of errors to allow before quitting\n";
      std::cout << "\t                       reading/validating the file.\n";
      std::cout << "\t                       -1 (default) indicates to not quit until\n";
      std::cout << "\t                       the entire file is read.\n";
      std::cout << "\t                       0 indicates not to read/validate anything\n";
      std::cout << "\t--printableErrors    : Maximum number of errors to print before\n";
      std::cout << "\t                       suppressing them (Defaults to 20).\n";
      std::cout << "\t                       Different than maxErrors since \n";
      std::cout << "\t                       printableErrors will continue reading and\n";
      std::cout << "\t                       validating the file until the end, but\n";
      std::cout << "\t                       just doesn't print the errors.\n";
      std::cout << "\t--ignoreErrors       : Ignore all errors (same as printableErrors = 0)\n";
      std::cout << "\t                       overwrites the printableErrors option.\n";
      std::cout << "\t--baseComposition    : Print the Base Composition Statistics.\n";
      std::cout << "\t--disableSeqIDCheck  : Disable the unique sequence identifier check.\n";
      std::cout << "\t                       Use this option to save memory since the sequence id\n";
      std::cout << "\t                       check uses a lot of memory.\n";
      std::cout << "\t--params             : Print the parameter settings.\n";
      std::cout << "\t--quiet              : Suppresses the display of errors and summary statistics.\n";
      std::cout << "\t                       Does not affect the printing of Base Composition Statistics.\n";

      std::cout << "\n  Optional Space Options for Raw Sequence (Last one specified is used):\n";
      std::cout << "\t--auto       : Determine baseSpace/colorSpace from the Raw Sequence in the file (Default).\n";
      std::cout << "\t--baseSpace  : ACTGN only\n";
      std::cout << "\t--colorSpace : 0123. only\n";
      std::cout << std::endl;

      std::cout << "  Usage:" << std::endl;
      std::cout << "\t./fastQValidator --file <fileName> [--minReadLen <minReadLen>] [--maxErrors <numErrors>] [--printableErrors <printableErrors>|--ignoreErrors] [--baseComposition] [--disableSeqIDCheck] [--quiet] [--baseSpace|--colorSpace|--auto] [--params]\n\n";
      std::cout << "  Examples:" << std::endl;
      std::cout << "\t../fastQValidator --file testFile.txt\n";
      std::cout << "\t../fastQValidator --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 100\n";
      std::cout << "\t./fastQValidator --file test/testFile.txt --minReadLen 10 --colorSpace --ignoreErrors\n";
      std::cout << std::endl;
      return (-1);
   }
   
   FastQFile validator(minReadLength, printableErrors);
   
   if(quiet)
   {
      validator.disableMessages();
   }

   if(disableSeqIDCheck)
   {
       validator.disableSeqIDCheck();
   }

   validator.setMaxErrors(maxErrors);

   FastQStatus::Status status = validator.validateFastQFile(filename, baseComposition, myBaseType);

   if(!quiet)
   {
      std::cout << "Returning: " << status << " : " << FastQStatus::getStatusString(status)
                << std::endl;
   }

   return(status);
}
