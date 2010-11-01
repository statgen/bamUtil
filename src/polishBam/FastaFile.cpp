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

#include "FastaFile.h"

FastaFile::FastaFile() : pcLine(NULL), nMaxLineLength(nDefaultMaxLineLength) {
}

FastaFile::~FastaFile() {
  close();
}

bool FastaFile::open(const char* fileName) {
  if ( sFileName.empty() && ( pcLine == NULL ) ) {
    sFileName = fileName;
    ifsIn.open(fileName, std::ios::in);
    if ( ifsIn.good() ) {
      nCurrentLine = 0;
      curSeqIndex = -1;
      curSeqLength = 0;
      allocateBuffer();
      return true;
    }
    else {
      return false;
    }
  }
  else {
    return false;
  }
}

void FastaFile::close() {
  //freeBuffer();
}

void FastaFile::allocateBuffer(int bufLength) {
  freeBuffer();
  pcLine = new char[bufLength];
  nMaxLineLength = bufLength;
  //std::cerr << "bufLength = " << bufLength << std::endl;
}

void FastaFile::freeBuffer() {
  if ( pcLine != NULL ) {
    delete[] pcLine;
  }
}

bool FastaFile::readLine() {
  int n;

  ifsIn.getline(pcLine, nMaxLineLength);
  if ( pcLine[0] == '>' ) {
    if ( curSeqIndex >= 0 ) {
      updateCurrentSequenceInfo();
    }
    std::vector<std::string> tokens;
    tokenizeString(pcLine, tokens);
    curSeqName = tokens[0].substr(1); // make SeqName without leading '>'
    MD5_Init(&curMD5Ctx);
    curSeqLength = 0;
    ++curSeqIndex;
  }
  else {
    n = strlen(pcLine);
    MD5_Update(&curMD5Ctx, pcLine, n);
    curSeqLength += n;
  }
  ++nCurrentLine;
  return ifsIn.good();
}

void FastaFile::readThru() {
  while(readLine());
  updateCurrentSequenceInfo();
}

void FastaFile::updateCurrentSequenceInfo() {
  unsigned char md5digest[1024];
  char md5string[1024];
  vsSequenceNames.push_back(curSeqName);
  vnSequenceLengths.push_back(curSeqLength);
  MD5_Final(md5digest, &curMD5Ctx);
  for(int i=0; i < 16; ++i) {
    sprintf(&md5string[i*2],"%02x",static_cast<int>(md5digest[i]));
  }
  vsMD5sums.push_back(md5string);
  //++curSeqIndex;
  //std::cerr << curSeqIndex << "\tSeqName=" << curSeqName << "\tSeqLength=" << curSeqLength << "\tMD5=" << md5string << std::endl;
}

uint32_t FastaFile::tokenizeString(const char* s, std::vector<std::string>& tokens) {
  if ( !tokens.empty() ) {
    tokens.clear();
  }

  // Assuming delimiters are space of tab, 
  // use stringstream to tokenize
  // In order to accept other delimiters, 
  // this routine needs to be rewritten
  std::stringstream ss (std::stringstream::in | std::stringstream::out);
  std::string tok;
  
  ss << s;
  
  while( ss >> tok ) {
    tokens.push_back(tok); // each token is considerd as strings
  }
  
  return(tokens.size()); // return the number of tokenized elements
}

/*
int main(int argc, char** argv) {
  FastaFile f;
  if ( f.open(argv[1]) ) {
    f.readThru();
    std::cerr << "CurrentLine = " << f.nCurrentLine << std::endl;
  }
  else {
    std::cerr << "Failed reading file" << std::endl;
  }
  return 0;
  }*/
