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

#ifndef __UNCOMPRESSEDFILETYPE_H__
#define __UNCOMPRESSEDFILETYPE_H__

#include <iostream>
#include <stdio.h>
#include "FileType.h"

class UncompressedFileType : public FileType
{
public:
    UncompressedFileType()
    {
        filePtr = NULL;
    }

    virtual ~UncompressedFileType()
    {
        if (filePtr != NULL)
        {
            close();
        }
    }

    UncompressedFileType(const char * filename, const char * mode);

    bool operator == (void * rhs)
    {
        // No two file pointers are the same, so if rhs is not NULL, then
        // the two pointers are different (false).
        if (rhs != NULL)
            return false;
        return (filePtr == rhs);
    }

    bool operator != (void * rhs)
    {
        // No two file pointers are the same, so if rhs is not NULL, then
        // the two pointers are different (true).
        if (rhs != NULL)
            return true;
        return (filePtr != rhs);
    }

    // Close the file.
    inline int close()
    {
        if((filePtr != stdout) && (filePtr != stdin))
        {
            int result = fclose(filePtr);
            filePtr = NULL;
            return result;
        }
        filePtr = NULL;
        return 0;
    }


    // Reset to the beginning of the file.
    inline void rewind()
    {
        // Just call rewind to move to the beginning of the file.
        ::rewind(filePtr);
    }

    // Check to see if we have reached the EOF.
    inline int eof()
    {
        //  check the file for eof.
        return feof(filePtr);
    }

    // Check to see if the file is open.
    virtual inline bool isOpen()
    {
        if (filePtr != NULL)
        {
            // filePtr is not null, so the file is open.
            return(true);
        }
        return(false);
    }

    // Write to the file
    inline unsigned int write(const void * buffer, unsigned int size)
    {
        return fwrite(buffer, 1, size, filePtr);
    }

    // Read into a buffer from the file.  Since the buffer is passed in and
    // this would bypass the fileBuffer used by this class, this method must
    // be protected.
    inline int read(void * buffer, unsigned int size)
    {
        return fread(buffer, 1, size, filePtr);
    }


    // Get current position in the file.
    // -1 return value indicates an error.
    virtual inline long int tell()
    {
        return ftell(filePtr);
    }


    // Seek to the specified offset from the origin.
    // origin can be any of the following:
    // Note: not all are valid for all filetypes.
    //   SEEK_SET - Beginning of file
    //   SEEK_CUR - Current position of the file pointer
    //   SEEK_END - End of file
    // Returns true on successful seek and false on a failed seek.
    virtual inline bool seek(long int offset, int origin)
    {
        long int returnVal = fseek(filePtr, offset, origin);
        // Check for success - 0 return value.
        if (returnVal == 0)
        {
            return true;
        }
        // Successful.
        return false;
    }


protected:
    // A FILE Pointer is used.
    FILE* filePtr;
};

#endif


