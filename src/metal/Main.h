/*
 *  Copyright (C) 2000-2010  Regents of the University of Michigan
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

#ifndef __METAL_H__
#define __METAL_H__

#include "StringArray.h"
#include "StringHash.h"
#include "MathVector.h"
#include "MathStats.h"
#include "IntArray.h"

class FileSummary
   {
   public:
      FileSummary(FileSummary * pointer = NULL)
         { next = pointer; }

      ~FileSummary()
         { if (next != NULL) delete next; }

      String filename;
      String header;
      String separators;

      double weight;
      double genomicControl;
      bool   useStrand;

      int markerColumn;
      int weightColumn;
      int pvalueColumn;
      int effectColumn;
      int firstColumn;
      int secondColumn;
      int stderrColumn;
      int inteffectColumn;
      int intstderrColumn;
      int intcovColumn;
      int freqColumn;
      int strandColumn;
      int minColumns;
      int expectedColumns;

      bool strictColumnCounting;
      bool logTransform;

      StringArray filterLabel;
      IntArray    filterColumn;
      IntArray    filterCondition;
      Vector      filterValue;
      StringArray filterAlternate;
      StringHash  filterSets;
      IntArray    filterCounts;

      int processedMarkers;

      FileSummary * next;
   };

// Crash and Control-C handlers
void UserBreak(int);
void OutOfMemory(int);
void SetupCrashHandlers();

// Basic setup
void ClearAll();

// Custom input filtering
#define LESS_THAN              1
#define LESS_THAN_OR_EQUAL     2
#define EQUAL_TO               3
#define GREATER_THAN_OR_EQUAL  4
#define GREATER_THAN           5
#define NOT_EQUAL_TO           6
#define STRING_MATCH           (EQUAL_TO | 0x80)
#define STRING_MISMATCH        (NOT_EQUAL_TO | 0x80)
#define STRING_IN_SET          (0x87)

void AddFilter(StringArray & filter);
int  TranslateCondition(String & condition);
void SetupFilters(StringArray & header);
bool ApplyFilter(StringArray & row);
void FilterSummary();
void ClearFilters();

// Processing of Input Files
void NumbersToLetters(String & al);
void FlipAllele(String & al);
bool FlipAlleles(String & al1, String & al2, double & effect, double & freq);
bool FlipIntAlleles(String & al1, String & al2, double & effect, double & inteffect, double & freq);
bool GuessSecondAllele(int marker, String & al1, String & al2);

// Workhorses
void Analyze(bool heterogeneity);
void ProcessFile(String & filename, FileSummary * history);
bool ReProcessFile(FileSummary * history);

// Help!
void ShowHelp(bool startup = false);

// General run control
void RunScript(FILE * file);


#endif

 
