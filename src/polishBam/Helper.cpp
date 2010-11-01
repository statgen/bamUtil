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

#include "Helper.h"
#include <string.h>

using namespace std;

void Logger::open(const char* filename, bool verbose) 
{
  b_verbose = verbose;
  if ( ( strlen(filename) > 0 ) && strcmp(filename,"__NONE__") != 0 ) {
    fp_log = fopen(filename, "w");
    if ( fp_log == NULL ) {
      fprintf(stderr,"ERROR: Cannot open the log file %s. Check if the directory exists and you have the permission to create a file", filename);
      abort();
    }
  }
  else {
    fp_log = NULL;
  }
  fp_err = stderr;
}

void Logger::write_log(const char* format, ... ) {
  va_list args;

  if ( fp_log != NULL ) {
    va_start (args, format);
    vfprintf(fp_log, format, args);
    va_end (args);
    fprintf(fp_log, "\n");
  }
  
  if ( b_verbose ) {
    va_start (args, format);
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
  }
}

void Logger::error(const char* format, ... ) {
  va_list args;

  if ( fp_log != NULL ) {
    va_start (args, format);
    fprintf(fp_log, "ERROR: ");
    vfprintf(fp_log, format, args);
    va_end (args);
    fprintf(fp_log, "\n");
  }

  va_start (args, format);
  fprintf(fp_err, "ERROR : ");
  vfprintf(fp_err, format, args);
  va_end (args);
  fprintf(fp_err, "\n");

  abort();
}

void Logger::warning(const char* format, ... ) {
  va_list args;

  if ( fp_log != NULL ) {
    va_start (args, format);
    fprintf(fp_log, "WARNING: ");
    vfprintf(fp_log, format, args);
    va_end (args);
    fprintf(fp_log, "\n");
  }

  va_start (args, format);
  fprintf(fp_err, "WARNING : ");
  vfprintf(fp_err, format, args);
  va_end (args);
  fprintf(fp_err, "\n");
}
