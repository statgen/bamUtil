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

#ifndef __INPUTFILE_H__
#define __INPUTFILE_H__

#ifdef  __gnu_linux__
#ifndef __ZLIB_AVAILABLE__
#define __ZLIB_AVAILABLE__
#endif
#endif

#include <stdio.h>
#include <iostream>
#include <cstring>

#include "FileType.h"

class InputFile
{
public:

    // DEFAULT means to use the default type to open a file for write.
    // The default type is equivalent to UNCOMPRESSED.
    enum ifileCompression {DEFAULT, UNCOMPRESSED, GZIP, BGZF};

    InputFile()
    {
        myFileTypePtr = NULL;
        myBufferIndex = 0;
        myCurrentBufferSize = 0;
        // Default to buffer.
        myAllocatedBufferSize = DEFAULT_BUFFER_SIZE;
        myFileBuffer = new char[myAllocatedBufferSize];
        myFileName.clear();
    }

    // Destructor
    ~InputFile();

    InputFile(const char * filename, const char * mode,
              InputFile::ifileCompression compressionMode = InputFile::DEFAULT);

    /// Buffering reads disables the tell call for bgzf files.
    /// Default buffer size is 1048576.
    /// Any previous values in the buffer will be deleted.
    /// Turn off read buffering by setting bufferSize = 1;
    inline void bufferReads(unsigned int bufferSize = DEFAULT_BUFFER_SIZE)
    {
        // If the buffer size is the same, do nothing.
        if(bufferSize == myAllocatedBufferSize)
        {
            return;
        }
        // Delete the previous buffer.
        if(myFileBuffer != NULL)
        {
            delete[] myFileBuffer;
        }
        myBufferIndex = 0;
        myCurrentBufferSize = 0;
        // The buffer size must be at least 1 so one character can be
        // read and ifgetc can just assume reading into the buffer.
        if(bufferSize < 1)
        {
            bufferSize = 1;
        }
        myFileBuffer = new char[bufferSize];
        myAllocatedBufferSize = bufferSize;

        if(myFileTypePtr != NULL)
        {
            if(bufferSize == 1)
            {
                myFileTypePtr->setBuffered(false);
            }
            else
            {
                myFileTypePtr->setBuffered(true);
            }
        }
    }


    /// Disable read buffering.
    inline void disableBuffering()
    {
        bufferReads(1);
        if(myFileTypePtr != NULL)
        {
            myFileTypePtr->setBuffered(false);
        }
    }

    
    // Close the file.
    inline int ifclose()
    {
        if (myFileTypePtr == NULL)
        {
            return EOF;
        }
        int result = myFileTypePtr->close();
        delete myFileTypePtr;
        myFileTypePtr = NULL;
        myFileName.clear();
        return result;
    }

    inline int ifread(void * buffer, unsigned int size)
    {
        // There are 2 cases:
        //  1) There are already size available bytes in buffer.
        //  2) There are not size bytes in buffer.

        // Determine the number of available bytes in the buffer.
        unsigned int availableBytes = myCurrentBufferSize - myBufferIndex;
        unsigned int returnSize = 0;

        // Case 1: There are already size available bytes in buffer.
        if (size <= availableBytes)
        {
            //   Just copy from the buffer, increment the index and return.
            memcpy(buffer, myFileBuffer+myBufferIndex, size);
            // Increment the buffer index.
            myBufferIndex += size;
            returnSize = size;
        }
        // Case 2: There are not size bytes in buffer.
        else
        {
            // Check to see if there are some bytes in the buffer.
            if (availableBytes > 0)
            {
                // Size > availableBytes > 0
                // Copy the available bytes into the buffer.
                memcpy(buffer, myFileBuffer+myBufferIndex, availableBytes);
            }
            unsigned int remainingSize = size - availableBytes;

            // Check if the remaining size is more or less than the
            // max buffer size.
            if(remainingSize < myAllocatedBufferSize)
            {
                // the remaining size is not the full buffer, but read
                //  a full buffer worth of data anyway.
                myCurrentBufferSize =
                    readFromFile(myFileBuffer, myAllocatedBufferSize);
                
                // Check to see how much was copied.
                unsigned int copySize = remainingSize;
                if(copySize > myCurrentBufferSize)
                {
                    copySize = myCurrentBufferSize;
                }

                // Now copy the rest of the bytes into the buffer.
                memcpy((char*)buffer+availableBytes, myFileBuffer, copySize);

                // set the buffer index to the location after what we read.
                myBufferIndex = copySize;
                
                returnSize = availableBytes + copySize;
            }
            else
            {
                // More remaining to be read than the max buffer size, so just
                // read directly into the output buffer.
                int readSize = readFromFile((char*)buffer + availableBytes,
                                            remainingSize);
                returnSize = readSize + availableBytes;
            }
        }
        return(returnSize);
    }


    // Get a character from the file.  Read a character from the internal
    // buffer, or if the end of the buffer has been reached, read from the
    // file into the buffer and return index 0.
    inline int ifgetc()
    {
        if (myBufferIndex >= myCurrentBufferSize)
        {
            // at the last index, read a new buffer.
            myCurrentBufferSize = readFromFile(myFileBuffer, myAllocatedBufferSize);
            myBufferIndex = 0;
        }
        // If the buffer index is still greater than or equal to the
        // myCurrentBufferSize, then we failed to read the file - return EOF.
        if (myBufferIndex >= myCurrentBufferSize)
        {
            return(EOF);
        }
        return(myFileBuffer[myBufferIndex++]);
    }

    // Reset to the beginning of the file.
    inline void ifrewind()
    {
        // Just set the myBufferIndex and the myCurrentBufferSize to 0 to simulate
        // clearing the buffer and call rewind to move to the beginning of the
        // file.
        if (myFileTypePtr == NULL)
        {
            // No pointer, so nothing to rewind.
            return;
        }
        myCurrentBufferSize = 0;
        myBufferIndex = 0;
        myFileTypePtr->rewind();
    }


    // Check to see if we have reached the EOF.
    inline int ifeof()
    {
        // Not EOF if we are not at the end of the buffer.
        if (myBufferIndex < myCurrentBufferSize)
        {
            // There are still available bytes in the buffer, so NOT EOF.
            return false;
        }
        else
        {
            if (myFileTypePtr == NULL)
            {
                // No myFileTypePtr, so not eof (return 0).
                return 0;
            }
            // exhausted our buffer, so check the file for eof.
            return myFileTypePtr->eof();
        }
    }

    // We do not buffer the write call, so just leave this as normal.
    inline unsigned int ifwrite(const void * buffer, unsigned int size)
    {
        if (myFileTypePtr == NULL)
        {
            // No myFileTypePtr, so return 0 - nothing written.
            return 0;
        }
        return myFileTypePtr->write(buffer, size);
    }

    // Returns whether or not the file was successfully opened.
    inline bool isOpen()
    {
        // It is open if the myFileTypePtr is set and says it is open.
        if ((myFileTypePtr != NULL) && myFileTypePtr->isOpen())
        {
            return true;
        }
        // File was not successfully opened.
        return false;
    }

    // Get current position in the file.
    // -1 return value indicates an error.
    inline long int iftell()
    {
        if (myFileTypePtr == NULL)
        {
            // No myFileTypePtr, so return false - could not seek.
            return -1;
        }
        return myFileTypePtr->tell();
    }


    // Seek to the specified offset from the origin.
    // origin can be any of the following:
    // Note: not all are valid for all filetypes.
    //   SEEK_SET - Beginning of file
    //   SEEK_CUR - Current position of the file pointer
    //   SEEK_END - End of file
    // Returns true on successful seek and false on a failed seek.
    inline bool ifseek(long int offset, int origin)
    {
        if (myFileTypePtr == NULL)
        {
            // No myFileTypePtr, so return false - could not seek.
            return false;
        }
        // Reset buffering since a seek is being done.
        myBufferIndex = 0;
        myCurrentBufferSize = 0;
        return myFileTypePtr->seek(offset, origin);
    }

    const char* getFileName() const
    {
        return(myFileName.c_str());
    }

protected:
    // Open a file. Called by the constructor.
    // Returns true if the file was successfully opened, false otherwise.
    bool openFile(const char * filename, const char * mode,
                  InputFile::ifileCompression compressionMode);

    // Read into a buffer from the file.  Since the buffer is passed in and
    // this would bypass the myFileBuffer used by this class, this method must
    // be protected.
    inline int readFromFile(void * buffer, unsigned int size)
    {
        // If no myFileTypePtr, return 0 - nothing read.
        if (myFileTypePtr == NULL)
        {
            return 0;
        }
        return myFileTypePtr->read(buffer, size);
    }

#ifdef __ZLIB_AVAILABLE__
    // Only necessary with zlib to determine what file type on a new
    // file.  Without zlib, there are only uncompressed files, so a special
    // method is not needed to determine the type of file to open.
    // Open a file.  This method will open a file with the specified name and
    // mode with the fileTypePtr associated with the specified compressionMode.
    void openFileUsingMode(const char* filename, const char* mode,
                           InputFile::ifileCompression compressionMode);
#endif

    // The size of the buffer used by this class.
    static const unsigned int DEFAULT_BUFFER_SIZE = 1048576;

    // Pointer to a class that interfaces with different file types.
    FileType* myFileTypePtr;

    unsigned int myAllocatedBufferSize;

    // Buffer used to do large reads rather than 1 by 1 character reads
    // from the file.  The class is then managed to iterate through the buffer.
    char* myFileBuffer;

    // Current index into the buffer.  Used to track where we are in reading the
    // file from the buffer.
    unsigned int myBufferIndex;

    // Current number of entries in the buffer.  Used to ensure that
    // if a read did not fill the buffer, we stop before hitting the
    // end of what was read.
    unsigned int myCurrentBufferSize;

    std::string myFileName;
};

typedef InputFile* IFILE;




inline IFILE ifopen(const char * filename, const char * mode,
                    InputFile::ifileCompression compressionMode = InputFile::DEFAULT)
{
    IFILE file = new InputFile(filename, mode, compressionMode);
    if (!file->isOpen())
    {

        // Not open, so delete the file, and return null.
        delete file;
        file = NULL;
    }
    return file;
}


inline int ifclose(IFILE file)
{
    int result = file->ifclose();
    delete file;
    file = NULL;
    return(result);
}

inline unsigned int ifread(IFILE file, void * buffer, unsigned int size)
{
    return(file->ifread(buffer, size));
}

inline int ifgetc(IFILE file)
{
    return(file->ifgetc());
}

inline void ifrewind(IFILE file)
{
    file->ifrewind();
}

inline int ifeof(IFILE file)
{
    return(file->ifeof());
}

inline unsigned int ifwrite(IFILE file, const void * buffer, unsigned int size)
{
    return(file->ifwrite(buffer, size));
}

// Get current position in the file.
// -1 return value indicates an error.
inline long int iftell(IFILE file)
{
    return (file->iftell());
}

// Seek to the specified offset from the origin.
// origin can be any of the following:
// Note: not all are valid for all filetypes.
//   SEEK_SET - Beginning of file
//   SEEK_CUR - Current position of the file pointer
//   SEEK_END - End of file
// Returns true on successful seek and false on a failed seek.
inline bool ifseek(IFILE file, long int offset, int origin)
{
    return (file->ifseek(offset, origin));
}

int ifprintf(IFILE output, char * format, ...);

inline IFILE operator >> (IFILE stream, std::string &str)
{
    str.clear();
    int ch;
    // not safe... newline handling?
    while ((ch = stream->ifgetc())!=EOF && (ch != '\n')) str.push_back(ch);
    return stream;
}

#endif

