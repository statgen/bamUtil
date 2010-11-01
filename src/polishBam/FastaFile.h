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

#ifndef __FASTA_FILE_H
#define __FASTA_FILE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <limits.h>
#include <stdint.h>
#include <openssl/md5.h>

class FastaFile {
public:
  std::string sFileName;
  std::vector<std::string> vsSequenceNames;
  std::vector<std::string> vsMD5sums;
  std::vector<uint32_t> vnSequenceLengths;
  std::ifstream ifsIn;
  char* pcLine;
  static const uint32_t nDefaultMaxLineLength = 100000;
  uint32_t nMaxLineLength;

  MD5_CTX curMD5Ctx;
  std::string curSeqName;
  uint32_t curSeqLength;
  int curSeqIndex;
  uint32_t nCurrentLine;

  FastaFile();
  ~FastaFile();
  bool open(const char* fileName);
  void close();
  void allocateBuffer(int bufLength = nDefaultMaxLineLength);
  void freeBuffer();
  bool readLine();
  void readThru();
  void updateCurrentSequenceInfo();
  static uint32_t tokenizeString(const char* s, std::vector<std::string>& tokens);
};

#endif // __FASTA_FILE_H
