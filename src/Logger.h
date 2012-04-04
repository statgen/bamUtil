/*
 *  Copyright (C) 2010-2012  Hyun Min Kang, 
 *                           Regents of the University of Michigan
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

#ifndef __FASTA_FILE_H__
#define __FASTA_FILE_H__


#include <iostream>

// Logger class for logging/error/warning
class Logger {
 protected:
  FILE* fp_log;
  FILE* fp_err;
  bool b_verbose;

  Logger() {} // default constructor prohibited
 public:
  static Logger* gLogger;
  Logger(const char* filename, bool verbose);
  void writeLog(const char* format, ...);
  void error(const char* format, ...);
  void warning(const char* format, ...);
};

#endif
