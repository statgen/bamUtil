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
#include "InputFileTest.h"
#include <assert.h>
#include <iostream>

#ifdef __ZLIB_AVAILABLE__
void testWrite();


int main(int argc, char ** argv)
{

   IFILE_Test myFile;

   myFile.test();

   testWrite();
}


const int IFILE_Test::TEST_FILE_SIZE = 37;
const int IFILE_Test::BGZF_TEST_FILE_SIZE = 93;
const std::string IFILE_Test::TEST_FILE_CONTENTS = "ABCDabcd1234\nEFGefg567\nhijklHIJKL8910";

void IFILE_Test::test()
{
   std::cout << "\nUncompressedFileType Tests:" << std::endl;
   test_readFromFile("txt");
   test_ifeof_ifrewind("txt");
   test_ifread_ifgetc("txt");
   test_ifclose("txt");
   test_ifseek("txt");
   test_noExistRead("txt");

   std::cout << "\nGzipFileType Tests:" << std::endl;
   test_readFromFile("gz");
   test_ifeof_ifrewind("gz");
   test_ifread_ifgetc("gz");
   test_ifclose("gz");
   test_ifseek("gz");
   test_noExistRead("gz");

   std::cout << "\nBgzfFileType Tests:" << std::endl;
   test_readFromFile("bam");
   test_ifeof_ifrewind("bam");
   test_ifread_ifgetc("bam");
   test_ifclose("bam");
   test_ifseek("bam");
   test_noExistRead("bam");

   std::cout << "\n .glf file Tests:" << std::endl;
   test_readFromFile("glf");
   test_ifeof_ifrewind("glf");
   test_ifread_ifgetc("glf");
   test_ifclose("glf");
   test_ifseek("glf");
   test_noExistRead("glf");
}


void IFILE_Test::test_readFromFile(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   // Track how many bytes are read by each call.
   int numBytesRead = 0;

   // Track the total number of the bytes that have been read from the file
   // at any given point.
   int totalBytesPreviouslyRead = 0;

   // Test readFromFile.
   numBytesRead = readFromFile(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should not have affected the internal buffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   // Should not be at eof
   assert(myFileTypePtr->eof() == false);
   assert(ifeof() == false);

   // Read again to verify that the next characters could be read.
   numBytesRead = readFromFile(myTestBuffer, 2);
   // Read 2 more characters from the test file.
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[4]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[5]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 2);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should not have affected the internal buffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   // Should not be at eof
   assert(myFileTypePtr->eof() == false);
   assert(ifeof() == false);
  
   // Read the rest of the file.
   // Determine expected results for reading the rest of the file by
   // taking the substring starting after what had been previously read.
   numBytesRead = readFromFile(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   // Read the rest of the file, so the number of bytes read is
   // what was left in the file.
   assert(numBytesRead == (TEST_FILE_SIZE - totalBytesPreviouslyRead));
   for(int i = 0; i < numBytesRead; i++)
   {
      assert(myTestBuffer[i] ==
	     TEST_FILE_CONTENTS[totalBytesPreviouslyRead+i]);
   }
   totalBytesPreviouslyRead += numBytesRead;
   // Eof - slightly varies based on the type of file on whether or not the
   // next read is the one the current read or the next one shows eof.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
      assert(myFileTypePtr->eof() == false);
      assert(ifeof() == false);      
   }
   else
   {
      assert(myFileTypePtr->eof() == true);
      assert(ifeof() == true);
   }

   // Try to read one more time, making sure it doesn't read anything.
    numBytesRead = readFromFile(myTestBuffer, MAX_TEST_BUFFER_SIZE);
    assert(numBytesRead == 0);
   // Should be at eof
   assert(myFileTypePtr->eof() == true);
   assert(ifeof() == true);

   ifclose();

   std::cout << "  Passed test_readFromFile" << std::endl;
}


void IFILE_Test::test_ifeof_ifrewind(const char* extension)
{
   // First open the test file.
   openFile(extension);
   
   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());
   
   // Not at eof - verify that it reports not eof.
   assert(ifeof() == false);
   
   // Track the total number of the bytes that have been read from the file
   // at any given point.
   int totalBytesPreviouslyRead = 0;

   //////////////////////////////////////////////////////////////
   // Test doing reads from file without IFILE internal buffering.
   // Read the entire file, then check for eof.
   int numBytesRead = readFromFile(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead < MAX_TEST_BUFFER_SIZE);

   // Eof - slightly varies based on the type of file on whether or not the
   // next read is the one the current read or the next one shows eof.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
      assert(myFileTypePtr->eof() == false);
      assert(ifeof() == false);      
   }
   else
   {
      assert(myFileTypePtr->eof() == true);
      assert(ifeof() == true);
   }
   
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 0);
   // Now it registers eof
   assert(ifeof() == true);

   // bgzf files use a specialized return value for iftell that
   // is not just straight file offset.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
       bool caught = false;
       try
       {
           assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
       }
       catch (std::exception& e)
       {
           caught = true;
           assert(strcmp(e.what(), "IFILE: CANNOT use buffered reads and tell for BGZF files") == 0);
       }
       assert(caught);
       disableBuffering();
       assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
   }
   else
   {
      assert(iftell() == TEST_FILE_SIZE);
   }

   // rewind the file and verify that it no longer registers eof.
   ifrewind();
   totalBytesPreviouslyRead = 0;
   // No longer at eof
   assert(ifeof() == false);
   // Verify position in file.
   assert(iftell() == 0);
   
   // Buffer reads - may have been disabled for iftell to work for bgzf.
   bufferReads();

   // Read a character from the file.
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;
   // Not at eof
   assert(ifeof() == false);

   ///////////////////////////////////
   // Test doing IFILE buffered reads.
   // Perform char read.
   char readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;
   // Not at eof
   assert(ifeof() == false);
   
   // Now read the rest.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == (TEST_FILE_SIZE - totalBytesPreviouslyRead));
   // Now that we have tested based on the previous total bytes read, 
   // increment the count.
   totalBytesPreviouslyRead += numBytesRead;
   // Registers eof.
   assert(ifeof() == true);

   // Read past eof.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Eof.
   assert(ifeof() == true);

   // bgzf files use a specialized return value for iftell that
   // is not just straight file offset.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
       bool caught = false;
       try
       {
           assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
       }
       catch (std::exception& e)
       {
           caught = true;
           assert(strcmp(e.what(), "IFILE: CANNOT use buffered reads and tell for BGZF files") == 0);
       }
       assert(caught);
       disableBuffering();
       assert(iftell() == (BGZF_TEST_FILE_SIZE << 16));
   }
   else
   {
      assert(iftell() == TEST_FILE_SIZE);
   }

   // Verify that after rewind, eof is no longer registered.
   ifrewind();
  // reset since we are back to the beginning of the file.
   totalBytesPreviouslyRead = 0;
   // No longer at eof
   assert(ifeof() == false);
   // Verify position in file.
   assert(iftell() == 0);

   // Verify properly works even if already at the beginning.   
   ifrewind();
  // reset since we are back to the beginning of the file.
   totalBytesPreviouslyRead = 0;
   // Not eof
   assert(ifeof() == false);
   // Verify position in file.
   assert(iftell() == 0);

   // Buffer reads - may have been disabled for iftell to work for bgzf.
   bufferReads();

   //////////////////////
   // Close the test file.
   ifclose();
   
   std::cout << "  Passed test_ifeof_ifrewind" << std::endl;
}


void IFILE_Test::test_ifread_ifgetc(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   int numBytesRead = 0;
   int totalBytesPreviouslyRead = 0;

   ////////////////////////////////////
   // Test reading entire file at once.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == TEST_FILE_SIZE);
   
   for(int i = 0; i < TEST_FILE_SIZE; i++)
   {
      assert(myTestBuffer[i] == TEST_FILE_CONTENTS[i]);
   }
   totalBytesPreviouslyRead += numBytesRead;
  
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == (unsigned int)TEST_FILE_SIZE);
   assert(myBufferIndex == (unsigned int)TEST_FILE_SIZE);
   
   // Eof - slightly varies based on the type of file on whether or not the
   // next read is the one the current read or the next one shows eof.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
      // does not yet register eof - but will on the next call.
      assert(myFileTypePtr->eof() == false);
      assert(ifeof() == false);      
   }
   else
   {
      assert(myFileTypePtr->eof() == true);
      assert(ifeof() == true);
   }

   // Try reading at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);
   

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////////
   // Test reading entire file using getc.
   // Loop through reading the file.
   char readChar;
   for(int index = 0; index < TEST_FILE_SIZE; index++)
   {
      // Read a character.
      readChar = ifgetc();
      assert(readChar == TEST_FILE_CONTENTS[index]);
      // Should affect the IFILE buffer
      assert(myCurrentBufferSize == (unsigned int)TEST_FILE_SIZE);
      assert(myBufferIndex == index+1);
   }
   
   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // Try again at eof.
   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////////////////
   // Test reading just the beginning of the file.
   numBytesRead = ifread(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == (unsigned int)TEST_FILE_SIZE);
   assert(myBufferIndex == 4);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading rest of file.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == TEST_FILE_SIZE - totalBytesPreviouslyRead);
   // Verify contents of what read.   
    for(int i = 0; i < numBytesRead; i++)
    {
       assert(myTestBuffer[i] == 
	      TEST_FILE_CONTENTS[i + totalBytesPreviouslyRead]);
    }
    totalBytesPreviouslyRead += numBytesRead;

   // Try at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////
   // Test reading just the beginning.
   numBytesRead = ifread(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == (unsigned int)TEST_FILE_SIZE);
   assert(myBufferIndex == 4);
   // Should not be at eof
   assert(ifeof() == false);

   // Test doing 2 getc.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   int bufferSize = TEST_FILE_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 5);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 6);
   totalBytesPreviouslyRead++;

   // Test reading rest of file.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == TEST_FILE_SIZE - totalBytesPreviouslyRead);
   // Verify contents of what read.   
    for(int i = 0; i < numBytesRead; i++)
    {
       assert(myTestBuffer[i] == 
	      TEST_FILE_CONTENTS[i + totalBytesPreviouslyRead]);
    }
    totalBytesPreviouslyRead += numBytesRead;

   // Try at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   //////////////////////////////////   
   // Start with 2 getc.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   bufferSize = TEST_FILE_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 1);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 2);
   totalBytesPreviouslyRead++;

   // Test reading part of the rest of the file.
   numBytesRead = ifread(myTestBuffer, 4);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myTestBuffer[1] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead + 1]);
   assert(myTestBuffer[2] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead + 2]);
   assert(myTestBuffer[3] == TEST_FILE_CONTENTS[totalBytesPreviouslyRead + 3]);
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == totalBytesPreviouslyRead);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading 2 char with getc.
   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   readChar = ifgetc();
   assert(readChar == TEST_FILE_CONTENTS[totalBytesPreviouslyRead]);
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   // Test reading rest of file.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == TEST_FILE_SIZE - totalBytesPreviouslyRead);
   // Verify contents of what read.   
    for(int i = 0; i < numBytesRead; i++)
    {
       assert(myTestBuffer[i] == 
	      TEST_FILE_CONTENTS[i + totalBytesPreviouslyRead]);
    }
    totalBytesPreviouslyRead += numBytesRead;
   assert(myBufferIndex == 0);
   assert(myCurrentBufferSize == 0);

   // Try at end of file twice.
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(myTestBuffer, MAX_TEST_BUFFER_SIZE);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   //////////////
   // Close the test file.
   ifclose();

   ////////////////////////////////////////////////////////////////////////
   // Repeat the test on a test file that is larger than the IFILE
   // buffer size.

   // First open the test file.
   openLargeFile(extension);

   // This file contains DEFAULT_BUFFER_SIZE of '0's followed by "12345"
   // The size of the file is DEFAULT_BUFFER_SIZE + 5.
   int largeTestFileSize = DEFAULT_BUFFER_SIZE + 5;
   char largeBuffer[largeTestFileSize + 5];

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   numBytesRead = 0;
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////
   // Test reading entire file at once.
   numBytesRead = ifread(largeBuffer, largeTestFileSize + 4);
   assert(numBytesRead == largeTestFileSize);
   
   // Validate all the 0s
   for(int i = 0; i < DEFAULT_BUFFER_SIZE; i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Now validate the "12345"
   assert(largeBuffer[DEFAULT_BUFFER_SIZE] == '1');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+1] == '2');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+2] == '3');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+3] == '4');
   assert(largeBuffer[DEFAULT_BUFFER_SIZE+4] == '5');

   totalBytesPreviouslyRead += numBytesRead;
  
   // Should affect the IFILE buffer - 0 because read
   // is bigger than the buffer, so just read directly
   // into the largeBuffer.
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   
   // Eof - slightly varies based on the type of file on whether or not the
   // next read is the one the current read or the next one shows eof.
   if((strcmp(extension, "bam") == 0) || (strcmp(extension, "glf") == 0))
   {
      // does not yet register eof - but will on the next call.
      assert(myFileTypePtr->eof() == false);
      assert(ifeof() == false);      
   }
   else
   {
      assert(myFileTypePtr->eof() == true);
      assert(ifeof() == true);
   }

   // Try reading at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);
   assert(ifeof() == true);
   

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////////
   // Test reading entire file using getc.
   // Loop through reading the file.
   // First loop through verifying the 0's
   for(unsigned int index = 0; index < DEFAULT_BUFFER_SIZE; index++)
   {
      // Read a character.
      readChar = ifgetc();
      assert(readChar == '0');
      // Should affect the IFILE buffer
      assert(myCurrentBufferSize == DEFAULT_BUFFER_SIZE);
      assert(myBufferIndex == index+1);
   }
   // Now read the 12345.
   readChar = ifgetc();
   assert(readChar == '1');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 1);
   readChar = ifgetc();
   assert(readChar == '2');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 2);
   readChar = ifgetc();
   assert(readChar == '3');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 3);
   readChar = ifgetc();
   assert(readChar == '4');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 4);
   readChar = ifgetc();
   assert(readChar == '5');
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 5);

   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // Try again at eof.
   // Now that we have read the file, try reading again at eof.
   readChar = ifgetc();
   assert(readChar == EOF);
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   ////////////////////////////////////////////////
   // Test reading just the beginning of the file.
   numBytesRead = ifread(largeBuffer, 4);
   assert(largeBuffer[0] == '0');
   assert(largeBuffer[1] == '0');
   assert(largeBuffer[2] == '0');
   assert(largeBuffer[3] == '0');
   assert(numBytesRead == 4);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == DEFAULT_BUFFER_SIZE);
   assert(myBufferIndex == 4);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading rest of file.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == largeTestFileSize - totalBytesPreviouslyRead);
   // Verify contents of what read.   First check the 0's
   for(int i = 0; i < (numBytesRead-5); i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Check the 12345
   assert(largeBuffer[numBytesRead - 5] == '1');
   assert(largeBuffer[numBytesRead - 5 + 1] == '2');
   assert(largeBuffer[numBytesRead - 5 + 2] == '3');
   assert(largeBuffer[numBytesRead - 5 + 3] == '4');
   assert(largeBuffer[numBytesRead - 5 + 4] == '5');
   totalBytesPreviouslyRead += numBytesRead;
   
   // Try at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer - size if 5, the ammount
   // that did not fit in the first read.
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 5);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 5);
   assert(ifeof() == true);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;

   //////////////////////////////////////
   // Test reading just the beginning.
   numBytesRead = ifread(largeBuffer, 2);
   assert(largeBuffer[0] == '0');
   assert(largeBuffer[1] == '0');
   assert(numBytesRead == 2);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == DEFAULT_BUFFER_SIZE);
   assert(myBufferIndex == 2);
   // Should not be at eof
   assert(ifeof() == false);

   // Test doing 2 getc.
   readChar = ifgetc();
   assert(readChar == '0');
   bufferSize = DEFAULT_BUFFER_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 3);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 4);
   totalBytesPreviouslyRead++;

   // Test reading rest of file.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == largeTestFileSize - totalBytesPreviouslyRead);
   // Verify contents of what read.   
   // All except the last 5 should be '0'
    for(int i = 0; i < numBytesRead - 5; i++)
    {
       assert(largeBuffer[i] == '0');
    }
    assert(largeBuffer[numBytesRead - 5] == '1');
    assert(largeBuffer[numBytesRead - 4] == '2');
    assert(largeBuffer[numBytesRead - 3] == '3');
    assert(largeBuffer[numBytesRead - 2] == '4');
    assert(largeBuffer[numBytesRead - 1] == '5');

    totalBytesPreviouslyRead += numBytesRead;

   // Try at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should not affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 5);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should not affect the IFILE buffer
   assert(myCurrentBufferSize == 5);
   assert(myBufferIndex == 5);
   assert(ifeof() == true);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   //////////////////////////////////   
   // Start with 2 getc.
   readChar = ifgetc();
   assert(readChar == '0');
   bufferSize = DEFAULT_BUFFER_SIZE;
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 1);
   totalBytesPreviouslyRead++;

   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == 2);
   totalBytesPreviouslyRead++;

   // Test reading part of the rest of the file.
   numBytesRead = ifread(myTestBuffer, 2);
   assert(myTestBuffer[0] == '0');
   assert(myTestBuffer[1] == '0');
   assert(numBytesRead == 2);
   totalBytesPreviouslyRead += numBytesRead;
   // This read should have affected the internal buffer.
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == totalBytesPreviouslyRead);
   // Should not be at eof
   assert(ifeof() == false);

   // Test reading 2 char with getc.
   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   readChar = ifgetc();
   assert(readChar == '0');
   assert(myCurrentBufferSize == bufferSize);
   totalBytesPreviouslyRead++;
   assert(myBufferIndex == totalBytesPreviouslyRead);

   // Test reading rest of file.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == largeTestFileSize - totalBytesPreviouslyRead);
   // Verify contents of what read.   
   for(int i = 0; i < numBytesRead - 5; i++)
   {
      assert(largeBuffer[i] == '0');
   }
   // Verify the 12345
   assert(largeBuffer[numBytesRead - 5] == '1');
   assert(largeBuffer[numBytesRead - 5 + 1] == '2');
   assert(largeBuffer[numBytesRead - 5 + 2] == '3');
   assert(largeBuffer[numBytesRead - 5 + 3] == '4');
   assert(largeBuffer[numBytesRead - 5 + 4] == '5');
   totalBytesPreviouslyRead += numBytesRead;
   bufferSize = 5;
   assert(myBufferIndex == bufferSize);
   assert(myCurrentBufferSize == bufferSize);

   // Try at end of file twice.
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should not affect the IFILE buffer
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == bufferSize);
   assert(ifeof() == true);

   // 2nd read attempt at eof.   
   numBytesRead = ifread(largeBuffer, largeTestFileSize);
   assert(numBytesRead == 0);
   // Should not affect the IFILE buffer
   assert(myCurrentBufferSize == bufferSize);
   assert(myBufferIndex == bufferSize);
   assert(ifeof() == true);
   
    // RESET
   ifrewind();
   totalBytesPreviouslyRead = 0;
   assert(myCurrentBufferSize == 0);
   assert(myBufferIndex == 0);

   ifclose();

   std::cout << "  Passed test_ifread_ifgetc" << std::endl;
}


// Test closing a file.
void IFILE_Test::test_ifclose(const char* extension)
{
   // First open the test file.
   openFile(extension);

   // Verify the file successfully opened.
   assert(myFileTypePtr != NULL);
   assert(isOpen());
   assert(myFileTypePtr->isOpen());

   ifclose();

   assert(myFileTypePtr == NULL);
   assert(isOpen() == false);

   std::cout << "  Passed test_ifclose" << std::endl;
}


void IFILE_Test::test_ifseek(const char* extension)
{
    disableBuffering();
   // First open the test file.
   openFile(extension);

   // Read a character from the file.
   int numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[0]);

   // Get the current position.
   long int currentPos = iftell();
   
   // Read the next character from the file.
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[1]);

   // Seek to just before the character that was just read and read again
   // Should be the same character.
   assert(ifseek(currentPos, SEEK_SET) == true);
   numBytesRead = readFromFile(myTestBuffer, 1);
   assert(numBytesRead == 1);
   assert(myTestBuffer[0] == TEST_FILE_CONTENTS[1]);
   
   ifclose();

   assert(myFileTypePtr == NULL);
   assert(isOpen() == false);

   // Buffer reads - may have been disabled for iftell to work for bgzf.
   bufferReads();

   std::cout << "  Passed test_ifseek" << std::endl;
}

void IFILE_Test::test_noExistRead(const char* extension)
{
   openNoExistFile(extension);

}


// Open a file for testing.
void IFILE_Test::openFile(const char* extension)
{
   std::string filename = "data/InputFileTest.";
   filename += extension;
   assert(InputFile::openFile(filename.c_str(), "rt", InputFile::DEFAULT) == true);
}

// Open a file for testing.
void IFILE_Test::openLargeFile(const char* extension)
{
   std::string filename = "data/InputFileTestLarge.";
   filename += extension;
   assert(InputFile::openFile(filename.data(), "rt", InputFile::DEFAULT) == true);
}


void IFILE_Test::openNoExistFile(const char* extension)
{
   std::string filename = "data/noExist.";
   filename += extension;
   assert(InputFile::openFile(filename.data(), "rt", InputFile::DEFAULT) == false);
}


void testWrite()
{
    std::string filenameNoExt = "results/InputFileTest.";
    std::string filename = filenameNoExt + "glf";
    
    IFILE filePtr = ifopen(filename.c_str(), "wt");
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr, 
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/uncompressedFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt", InputFile::UNCOMPRESSED);
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/bgzfFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt", InputFile::BGZF);
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/gzipFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt", InputFile::GZIP);
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           ==IFILE_Test:: TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/defaultFile.glf";
    
    filePtr = ifopen(filename.c_str(), "wt");
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    filename = "results/defaultFile.gz";
    
    filePtr = ifopen(filename.c_str(), "wt");
    assert(filePtr != NULL);
    
    assert(ifwrite(filePtr,
                   IFILE_Test::TEST_FILE_CONTENTS.c_str(), 
                   IFILE_Test::TEST_FILE_CONTENTS.length()) 
           == IFILE_Test::TEST_FILE_CONTENTS.length());
    
    assert(ifclose(filePtr) == 0);

    // TODO - automatically verify that the files were written in the
    // correct format - rather than hand checking.
}


#endif
