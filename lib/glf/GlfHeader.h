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

#ifndef __GLF_HEADER_H__
#define __GLF_HEADER_H__

#include <stdint.h>

#include "InputFile.h" 
#include "CharBuffer.h"

class GlfHeader
{
public:
    GlfHeader();
    ~GlfHeader();

    // Copy Constructor   
    GlfHeader(const GlfHeader& header);

    // Overload operator = to copy the passed in header into this header.
    GlfHeader & operator = (const GlfHeader& header);

    // Overload operator = to copy the passed in header into this header.
    bool copy(const GlfHeader& header);

    void resetHeader();
   
    // Read the header from the specified file.  Assumes the file is in
    // the correct position for reading the header.
    bool read(IFILE filePtr);

    // Write the header to the specified file.
    bool write(IFILE filePtr) const;

    // Set the passed in string to the text string stored in this header.
    bool getHeaderTextString(std::string& text);

    // Set the header to the passed in string.
    bool setHeaderTextString(const std::string& text);

private:
    //    std::string myText;

    int32_t myTextLen;
    CharBuffer myText;

    static const std::string GLF_MAGIC;
    static const int GLF_MAGIC_LEN;
};

#endif

