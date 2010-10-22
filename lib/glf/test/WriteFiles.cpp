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

#include "GlfFile.h"
#include "WriteFiles.h"
#include "Validate.h"

#include <assert.h>
#include <stdio.h>

void testWrite()
{
    GlfFile glfOut;
    
    assert(glfOut.openForWrite("results/MyTestOut.glf"));

    // Create a glf header.
    GlfHeader glfHeader;
   
    std::string headerString = "";
    std::string expectedHeaderString = "This is my test header.";

    // Write the header.
    assert(glfHeader.getHeaderTextString(headerString));
    assert(headerString == "");
    assert(glfHeader.setHeaderTextString(expectedHeaderString));
    assert(glfHeader.getHeaderTextString(headerString));
    assert(headerString == expectedHeaderString);
    assert(glfOut.writeHeader(glfHeader));

    ////////////////////////////////
    // Write the reference section.
    GlfRefSection glfSection;
    std::string refNameString = "";
    std::string expectedRefNameString = "This is my 1st Ref Name";
    uint32_t expectedRefLen = 100;
    // Check the default settings (no data has been set yet).
    assert(glfSection.getName(refNameString));
    assert(refNameString == "");
    assert(glfSection.getRefLen() == 0);

    // Set the reference name.
    assert(glfSection.setName(expectedRefNameString));
    // Check properly set.
    assert(glfSection.getName(refNameString));
    assert(refNameString == expectedRefNameString);
    assert(glfSection.getRefLen() == 0);

    // Set the reference sequence length.
    assert(glfSection.setRefLen(expectedRefLen));
    // Check properly set.
    assert(glfSection.getRefLen() == expectedRefLen);
    assert(glfSection.getName(refNameString));
    assert(refNameString == expectedRefNameString);

    // Write the reference section
    assert(glfOut.writeRefSection(glfSection));

    //////////////////////////////////////////////
    // Write a record of type 1.
    GlfRecord record;
    uint8_t expectedRecType1 = 1;
    uint8_t expectedRefBase1 = 2;
    uint32_t expectedOffset1 = 500;
    uint32_t expectedMinLk1 = 55;
    uint32_t expectedReadDepth1 = 31;
    uint8_t expectedRmsMapQ1 = 25;
    
    assert(record.setRecordType(expectedRecType1));
    assert(record.setRefBaseInt(expectedRefBase1));
    assert(record.setOffset(expectedOffset1));
    assert(record.setMinLk(expectedMinLk1));
    assert(record.setReadDepth(expectedReadDepth1));
    assert(record.setRmsMapQ(expectedRmsMapQ1));
    assert(glfOut.writeRecord(record));

    // Verify the settings of record 1.
    assert(record.getRecordType() == expectedRecType1);
    assert(record.getRefBase() == expectedRefBase1);
    assert(record.getOffset() == expectedOffset1);
    assert(record.getMinLk() == expectedMinLk1);
    assert(record.getReadDepth() == expectedReadDepth1);
    assert(record.getRmsMapQ() == expectedRmsMapQ1);

    //////////////////////////////////////////////
    // Write a record of type 2.
    record.reset();
    uint8_t expectedRecType2 = 2;
    uint8_t expectedRefBase2 = 1;
    uint32_t expectedOffset2 = 6;
    uint32_t expectedMinLk2 = 44;
    uint32_t expectedReadDepth2 = 66;
    uint8_t expectedRmsMapQ2 = 32;
    uint8_t expectedLkHom1 = 98;
    uint8_t expectedLkHom2 = 86;
    uint8_t expectedLkHet = 73;
    uint16_t expectedIndelLen1 = 0;
    uint16_t expectedIndelLen2 = 0;
    std::string expectedIndelSeq1 = "";
    std::string expectedIndelSeq2 = "";
    std::string indelSeq = "";

    assert(record.setRecordType(expectedRecType2));
    assert(record.setRefBaseInt(expectedRefBase2));
    assert(record.setOffset(expectedOffset2));
    assert(record.setMinLk(expectedMinLk2));
    assert(record.setReadDepth(expectedReadDepth2));
    assert(record.setRmsMapQ(expectedRmsMapQ2));
    assert(record.setLkHom1(expectedLkHom1));
    assert(record.setLkHom2(expectedLkHom2));
    assert(record.setLkHet(expectedLkHet));
    assert(glfOut.writeRecord(record));

    // Verify the settings of record 2.
    assert(record.getRecordType() == expectedRecType2);
    assert(record.getRefBase() == expectedRefBase2);
    assert(record.getOffset() == expectedOffset2);
    assert(record.getMinLk() == expectedMinLk2);
    assert(record.getReadDepth() == expectedReadDepth2);
    assert(record.getRmsMapQ() == expectedRmsMapQ2);
    assert(record.getLkHom1() == expectedLkHom1);
    assert(record.getLkHom2() == expectedLkHom2);
    assert(record.getLkHet() == expectedLkHet);
    assert(record.getIndel1(indelSeq) == expectedIndelLen1);
    assert(indelSeq == expectedIndelSeq1);
    assert(record.getIndel2(indelSeq) == expectedIndelLen2);
    assert(indelSeq == expectedIndelSeq2);

    //////////////////////////////////////////////
    // Write a record of type 0.
    record.reset();
    assert(glfOut.writeRecord(record));

    // Verify the settings of the types.
    assert(record.getRecordType() == 0);
    assert(record.getRefBase() == 0);

    // Close the file.
    glfOut.close();

    // TODO test exceptions when wrong accessors are used.

    // Validate the just written file.
    GlfFile glfIn;
    assert(glfIn.openForRead("results/MyTestOut.glf"));

    // Check the header string.
    assert(glfIn.readHeader(glfHeader));
    assert(glfHeader.getHeaderTextString(headerString));
    assert(headerString == expectedHeaderString);

    // Check the reference section.
    assert(glfIn.getNextRefSection(glfSection));
    assert(glfSection.getName(refNameString));
    assert(refNameString == expectedRefNameString);
    assert(glfSection.getRefLen() == expectedRefLen);

    // Check the record of type 1.
    assert(glfIn.getNextRecord(record));
    assert(record.getRecordType() == expectedRecType1);
    assert(record.getRefBase() == expectedRefBase1);
    assert(record.getOffset() == expectedOffset1);
    assert(record.getMinLk() == expectedMinLk1);
    assert(record.getReadDepth() == expectedReadDepth1);
    assert(record.getRmsMapQ() == expectedRmsMapQ1);

    //Check the record of type 2.
    assert(glfIn.getNextRecord(record));
    assert(record.getRecordType() == expectedRecType2);
    assert(record.getRefBase() == expectedRefBase2);
    assert(record.getOffset() == expectedOffset2);
    assert(record.getMinLk() == expectedMinLk2);
    assert(record.getReadDepth() == expectedReadDepth2);
    assert(record.getRmsMapQ() == expectedRmsMapQ2);
    assert(record.getLkHom1() == expectedLkHom1);
    assert(record.getLkHom2() == expectedLkHom2);
    assert(record.getLkHet() == expectedLkHet);
    assert(record.getIndel1(indelSeq) == expectedIndelLen1);
    assert(indelSeq == expectedIndelSeq1);
    assert(record.getIndel2(indelSeq) == expectedIndelLen2);
    assert(indelSeq == expectedIndelSeq2);

    // Check the record of type 0.
    // False, since there are no more records in this section.
    assert(glfIn.getNextRecord(record) == false);
    assert(record.getRecordType() == 0);
    assert(record.getRefBase() == 0);
}

