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

#include "InputFile.h"
#include "StringBasics.h"
#include "GzipHeader.h"
#include "BgzfFileType.h"
#include "GzipFileType.h"
#include "UncompressedFileType.h"

#include <stdarg.h>

InputFile::InputFile(const char * filename, const char * mode,
                     InputFile::ifileCompression compressionMode)
{
    myFileTypePtr = NULL;
    myBufferIndex = 0;
    myCurrentBufferSize = 0;
    myAllocatedBufferSize = DEFAULT_BUFFER_SIZE;
    myFileBuffer = new char[myAllocatedBufferSize];
    myFileName.clear();

    openFile(filename, mode, compressionMode);
}


#ifdef __ZLIB_AVAILABLE__

// Open a file. Called by the constructor.
// Returns true if the file was successfully opened, false otherwise.
bool InputFile::openFile(const char * filename, const char * mode,
                         InputFile::ifileCompression compressionMode)
{
    // If a file is for write, just open a new file.
    if (mode[0] == 'w' || mode[0] == 'W')
    {
        openFileUsingMode(filename, mode, compressionMode);
    }
    else
    {
        // Check if reading from stdin.
        if(strcmp(filename, "-") == 0)
        {
            // Reading from stdin, open it based on the 
            // compression mode.
            openFileUsingMode(filename, mode, compressionMode);
        }
        else
        {
            // Not from stdin, so determine the file type.

            // Open the file read only to determine file type.
            FILE* filePtr = fopen(filename, "r");
            
            // If the file could not be opened, either create a new one or
            // return failure.
            if (filePtr == NULL)
            {
                // If the mode is for read, then the file must exist, otherwise,
                // create a new file.
                if (mode[0] == 'r' || mode[0] == 'R')
                {
                    // File must exist.
                    if (myFileTypePtr != NULL)
                    {
                        delete myFileTypePtr;
                        myFileTypePtr = NULL;
                    }
                    // Return false, was not opened.
                    return false;
                }
                else
                {
                    openFileUsingMode(filename, mode, compressionMode);
                }
            }
            else
            {
                // File was successfully opened, so try to determine the
                // filetype from the file.
                // Read the file to see if it a gzip file.
                GzipHeader gzipHeader;
                bool isGzip = gzipHeader.readHeader(filePtr);
                
                // The file header has been read, so close the file, so it can
                // be re-opened as the correct type.
                fclose(filePtr);

                if (isGzip)
                {
                    // This file is a gzip file.
                    // Check to see if it is BGZF Compression.
                    if (gzipHeader.isBgzfFile())
                    {
                        // This file has BGZF Compression, so set the file
                        // pointer.
                        myFileTypePtr = new BgzfFileType(filename, mode);
                    }
                    else
                    {
                        // Not BGZF, just a normal gzip.
                        myFileTypePtr = new GzipFileType(filename, mode);
                   }
                }
                else
                {
                    // The file is a uncompressed, uncompressed file,
                    // so set the myFileTypePtr accordingly.
                    myFileTypePtr = new UncompressedFileType(filename, mode);
                }
            }
        }
    }
    if(myFileTypePtr == NULL)
    {
        return(false);
    }
    if (!myFileTypePtr->isOpen())
    {
        // The file was not opened, so delete the pointer and set to null.
        delete myFileTypePtr;
        myFileTypePtr = NULL;
        return false;
    }

    if(myAllocatedBufferSize == 1)
    {
        myFileTypePtr->setBuffered(false);
    }
    else
    {
        myFileTypePtr->setBuffered(true);
    }
    myFileName = filename;
    return true;
}


// Open a file.  This method will open a file with the specified name and
// mode with the fileTypePtr associated with the specified compressionMode.
void InputFile::openFileUsingMode(const char * filename, const char * mode,
                                  ifileCompression compressionMode)
{
    switch (compressionMode)
    {
        case GZIP:
            // Gzipped.
            myFileTypePtr = new GzipFileType(filename, mode);
            break;
        case BGZF:
            // BGZF compression.
            myFileTypePtr = new BgzfFileType(filename, mode);
            break;
        case UNCOMPRESSED:
            myFileTypePtr = new UncompressedFileType(filename, mode);
            break;
        case InputFile::DEFAULT:
        default:
            // Check the extension. If it is ".gz", treat as gzip.
            // otherwise treat it as UNCOMPRESSED.
            int lastchar = 0;
            while (filename[lastchar] != 0) lastchar++;
            if ((lastchar >= 3 &&
                    filename[lastchar - 3] == '.' &&
                    filename[lastchar - 2] == 'g' &&
                    filename[lastchar - 1] == 'z'))
            {
                // .gz files files should be gzipped.
                myFileTypePtr = new GzipFileType(filename, mode);
            }
            else
            {
                // Create an uncompressed file.
                myFileTypePtr = new UncompressedFileType(filename, mode);
            }
            break;
    }

    if(myFileTypePtr == NULL)
    {
        return;
    }
    if(myAllocatedBufferSize == 1)
    {
        myFileTypePtr->setBuffered(false);
    }
    else
    {
        myFileTypePtr->setBuffered(true);
    }
}

#else

// No zlib, so just treat all files as std files.
// Open a file. Called by the constructor.
// Returns true if the file was successfully opened, false otherwise.
bool InputFile::openFile(const char * filename, const char * mode)
{
    //  No zlib, so it is a uncompressed, uncompressed file.
    myFileTypePtr = new UncompressedFileType(filename, mode);

    if(myFileTypePtr == NULL)
    {
        return(false);
    }
    if (!myFileTypePtr->isOpen())
    {
        // The file was not opened, so delete the pointer and set to null.
        delete myFileTypePtr;
        myFileTypePtr = NULL;
        return false;
    }
    if(myAllocatedBufferSize == 1)
    {
        myFileTypePtr->setBuffered(false);
    }
    else
    {
        myFileTypePtr->setBuffered(true);
    }
    myFileName = filename;
    return true;
}

#endif


InputFile::~InputFile()
{
    delete myFileTypePtr;
    myFileTypePtr = NULL;

    if(myFileBuffer != NULL)
    {
        delete[] myFileBuffer;
        myFileBuffer = NULL;
    }
}


int ifprintf(IFILE output, char * format, ...)
{
    String buffer;

    va_list  ap;
    va_start(ap, format);

    buffer.vprintf(format, ap);

    va_end(ap);

    return ::ifwrite(output, (const char *) buffer, buffer.Length());
}


