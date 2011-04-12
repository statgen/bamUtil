#ifndef __CSG_LOGGER_H__
#define __CSG_LOGGER_H__

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "InputFile.h"
#include "StringBasics.h"

class Logger {
 protected:
  IFILE fp_log;
  FILE* fp_err;
  bool b_verbose;
  
  Logger() {} // default constructor prohibited
 public:
  static Logger* gLogger;

  Logger(const char* filename, bool verbose) {
    b_verbose = verbose;
    if ( strlen(filename) > 0 ) {
      //strcmp(filename,"__NONE__") != 0 ) {
      fp_log = ifopen(filename, "wb");
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

  void writeLog(const char* format, ...) {
    va_list args;
    if ( fp_log != NULL ) {
      String buffer;
      va_start (args, format);
      buffer.vprintf(format, args);
      va_end (args);
      ::ifwrite(fp_log, (const char*) buffer, buffer.Length());
      //fflush(fp_log);
    }
    
    if ( b_verbose ) {
      va_start (args, format);
      vfprintf(fp_err, format, args);
      va_end (args);
      fprintf(fp_err, "\n");
      fflush(fp_err);
    }
  }

  void error(const char* format, ...) {
    va_list args;
    if ( fp_log != NULL ) {
      String buffer;
      va_start (args, format);
      buffer.printf("ERROR: ");
      buffer.vprintf(format, args);
      va_end (args);
      ::ifwrite(fp_log, (const char*) buffer, buffer.Length());
      //fflush(fp_log);
    }
    
    va_start (args, format);
    fprintf(fp_err, "ERROR : ");
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
    fflush(fp_err);

    abort();    
  }

  void warning(const char* format, ...) {
    va_list args;
    if ( fp_log != NULL ) {
      String buffer;
      va_start (args, format);
      buffer.printf("WARNING: ");
      buffer.vprintf(format, args);
      va_end (args);
      ::ifwrite(fp_log, (const char*) buffer, buffer.Length());
      //fflush(fp_log);
    }
    
    va_start (args, format);
    fprintf(fp_err, "WARNING : ");
    vfprintf(fp_err, format, args);
    va_end (args);
    fprintf(fp_err, "\n");
    fflush(fp_err);
  }
};

#endif // __LOGGER_H__
