/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#ifndef __VCF_FILE_H__
#define __VCF_FILE_H__


#include "InputFile.h"

class VcfFile
{
public:
    struct VcfDataLine
    {
        std::string wholeLine; // includes '\n'
        std::string chromosome;
        int position;
    };
    
    VcfFile();
    ~VcfFile();
    
    bool openForRead(const char* filename);
    bool openForWrite(const char* filename);
    void close();
    bool isEof();

    // Read the file until a data line is found and set nextLine.
    bool getNextDataLine(VcfDataLine& nextLine);

    // Returns true if successfully wrote or if a null is passed in.
    bool writeLine(VcfDataLine* line);

    // Returns true if successfully wrote or if a null is passed in.
    bool writeLine(std::string& line);

private:

    IFILE myFilePtr;
};
#endif
