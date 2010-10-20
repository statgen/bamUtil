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

#include "UncompressedFileType.h"
#include <string.h>

UncompressedFileType::UncompressedFileType(const char * filename,
                                           const char * mode)
{
    // If the file is for write and is '-', then write to stdout.
    if(((mode[0] == 'w') || (mode[0] == 'W')) && 
       (strcmp(filename, "-") == 0))
    {
        // Write to stdout.
        filePtr = stdout;
    }
    else if(((mode[0] == 'r') || (mode[0] == 'R')) && 
       (strcmp(filename, "-") == 0))
    {
        // read from stdin
        filePtr = stdin;
    }
    else
    {
        // Open the file.
        filePtr = fopen(filename, mode);
    }
};


