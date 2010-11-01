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

#ifndef __HELPER_H__
#define __HELPER_H__

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

//using namespace std;

class Logger {
 protected:
  FILE* fp_log;
  FILE* fp_err;
  bool b_verbose;
  
 public:
 Logger(bool verbose=false) : fp_log(NULL), fp_err(stderr), b_verbose(verbose) {} // default constructor prohibited
  void open(const char* file, bool verbose);
    //Logger(const char* filename, bool verbose);
  void write_log(const char* format, ...);
  void error(const char* format, ...);
  void warning(const char* format, ...);
  //void error(const char* s);
  //void warning(const char* s);
  //void log(const char* s);
};

#endif // __HELPER_H__
