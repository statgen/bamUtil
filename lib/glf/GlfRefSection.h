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

#ifndef __GLF_REFSECTION_H__
#define __GLF_REFSECTION_H__

#include <stdint.h>

#include "InputFile.h" 
#include "CharBuffer.h"

class GlfRefSection
{
public:
    GlfRefSection();
    ~GlfRefSection();

    // Copy Constructor   
    GlfRefSection(const GlfRefSection& refSection);

    // Overload operator = to copy the passed in refSection into this refSection.
    GlfRefSection & operator = (const GlfRefSection& refSection);

    // Overload operator = to copy the passed in refSection into this refSection.
    bool copy(const GlfRefSection& refSection);

    void resetRefSection();
   
    // Read the refSection from the specified file.  Assumes the file is in
    // the correct position for reading the refSection.
    bool read(IFILE filePtr);

    // Write the refSection to the specified file.
    bool write(IFILE filePtr) const;

    /////////////
    // Accessors.

    //Get
    bool getName(std::string& name) const;
    uint32_t getRefLen() const;

    // Set
    bool setName(const std::string& name);
    bool setRefLen(uint32_t refLen);

    void print() const;

private:

    CharBuffer myRefName;
    uint32_t myRefLen;
};

#endif

