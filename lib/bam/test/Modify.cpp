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

#include "SamFile.h"
#include "Modify.h"

void testModify()
{
    modify modTest;
    modTest.testModify("testFiles/testSam.sam");
    modTest.testModify("testFiles/testBam.bam");
}


void modify::testModify(const char* filename)
{
    myFilename = filename;

    modifyPosition();
    modifyCigar();
}

void modify::modifyPosition()
{
    openAndRead1Rec();
   
    // Verify the initial bin.
    assert(samRecord.getBin() == 4681);

    // Change the position and verify that the bin is updated.
    assert(samRecord.set0BasedPosition(33768));

    // Verify the bin was updated.
    assert(samRecord.getBin() == 4683);
    assert(samRecord.get0BasedPosition() == 33768);
}


void modify::modifyCigar()
{
    openAndRead1Rec();
   
    // Verify the initial bin.
    assert(samRecord.getBin() == 4681);

    // Change the Cigar such that it modifies the bin.
    assert(samRecord.setCigar("33768M"));

    // Verify the bin was updated.
    assert(samRecord.getBin() == 585);
}


void modify::openAndRead1Rec()
{
    // Open the file for reading.   
    assert(samIn.OpenForRead(myFilename));

    // Read the sam header.
    assert(samIn.ReadHeader(samHeader));
   
    // Read the first record.   
    assert(samIn.ReadRecord(samHeader, samRecord));
}
