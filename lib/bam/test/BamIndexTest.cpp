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

#include "BamIndex.h"
#include "Validate.h"

#include <assert.h>

void testBamIndex()
{
    // Create a bam index.
    BamIndex bamIndex;
    bamIndex.readIndex("testFiles/sortedBam.bam.bai");

    // Get the chunks for reference id 1.
    Chunk testChunk;
    SortedChunkList chunkList;
    assert(bamIndex.getChunksForRegion(1, -1, -1, chunkList) == true);
    assert(!chunkList.empty());
    testChunk = chunkList.pop();
    assert(chunkList.empty());
    assert(testChunk.chunk_beg == 0x4e7);
    assert(testChunk.chunk_end == 0x599);

    // Get the chunks for reference id 0.
    assert(bamIndex.getChunksForRegion(0, -1, -1, chunkList) == true);
    assert(!chunkList.empty());
    testChunk = chunkList.pop();
    assert(chunkList.empty());
    assert(testChunk.chunk_beg == 0x360);
    assert(testChunk.chunk_end == 0x4e7);


    // Get the chunks for reference id 2.
    assert(bamIndex.getChunksForRegion(2, -1, -1, chunkList) == true);
    assert(!chunkList.empty());
    testChunk = chunkList.pop();
    assert(chunkList.empty());
    assert(testChunk.chunk_beg == 0x599);
    assert(testChunk.chunk_end == 0x5ea);

    // Get the chunks for reference id 3.
    // There isn't one for this ref id, but still successfully read the file,
    // so it should return true, but the list should be empty.
    assert(bamIndex.getChunksForRegion(3, -1, -1, chunkList) == true);
    assert(chunkList.empty());

    // Test reading an indexed bam file.
    SamFile inFile;
    assert(inFile.OpenForRead("testFiles/sortedBam.bam"));
    inFile.setSortedValidation(SamFile::COORDINATE);
    assert(inFile.ReadBamIndex("testFiles/sortedBam.bam.bai"));
    SamFileHeader samHeader;
    assert(inFile.ReadHeader(samHeader));
    SamRecord samRecord;

    // Section -1 = Ref *: 2 records (8 & 10 from testSam.sam that is reflected
    // in the validation.
    assert(inFile.SetReadSection(-1));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead8(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead10(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 2 = Ref 3: 1 records (9 from testSam.sam that is reflected
    // in the validation.
    assert(inFile.SetReadSection(2));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead9(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 0 = Ref 1: 5 records (3, 4, 1, 2, & 6 from testSam.sam that is
    // reflected in the validation.
    assert(inFile.SetReadSection(0));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead3(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead4(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead6(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 1 = Ref 2: 2 records (5 & 7 from testSam.sam that is reflected
    // in the validation.
    assert(inFile.SetReadSection(1));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead5(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead7(samRecord);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);

    // Section 3 to 22 (ref 4 - 23): 0 records.
    for(int i = 3; i < 23; i++)
    {
        assert(inFile.SetReadSection(i));
        assert(inFile.ReadRecord(samHeader, samRecord) == false);
    }


    // Set the read section.
    assert(inFile.SetReadSection("1", 1010, 1012));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 2);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 0);
    assert(samRecord.getNumOverlaps(1010, 1012) == 0);
    assert(samRecord.getNumOverlaps(1010, 1020) == 0);
    assert(samRecord.getNumOverlaps(1010, 1011) == 0);
    assert(samRecord.getNumOverlaps(1011, 1012) == 0);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
    assert(inFile.SetReadSection("1", 1010, 1020));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 5);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 0);
    assert(samRecord.getNumOverlaps(1010, 1012) == 0);
    assert(samRecord.getNumOverlaps(1010, 1020) == 0);
    assert(samRecord.getNumOverlaps(1010, 1011) == 0);
    assert(samRecord.getNumOverlaps(1011, 1012) == 0);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
    assert(inFile.SetReadSection("1", 1010, 1011));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 1);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
    assert(inFile.SetReadSection("1", 1011, 1012));
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead1(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 1);
    assert(samRecord.getNumOverlaps(1010, 1012) == 2);
    assert(samRecord.getNumOverlaps(1010, 1020) == 5);
    assert(samRecord.getNumOverlaps(1010, 1011) == 1);
    assert(samRecord.getNumOverlaps(1011, 1012) == 1);
    assert(inFile.ReadRecord(samHeader, samRecord));
    validateRead2(samRecord);
    assert(inFile.GetNumOverlaps(samRecord) == 0);
    assert(samRecord.getNumOverlaps(1010, 1012) == 0);
    assert(samRecord.getNumOverlaps(1010, 1020) == 0);
    assert(samRecord.getNumOverlaps(1010, 1011) == 0);
    assert(samRecord.getNumOverlaps(1011, 1012) == 0);
    assert(inFile.ReadRecord(samHeader, samRecord) == false);
           
}
