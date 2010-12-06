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

#ifndef __BAM_INTERFACE_H__
#define __BAM_INTERFACE_H__

#include "GenericSamInterface.h"

class BamInterface : public GenericSamInterface
{
public:
    BamInterface();
    ~BamInterface();
   
    // Reads the header section from the specified BAM file and stores it in
    // the passed in header.
    virtual SamStatus::Status readHeader(IFILE filePtr, SamFileHeader& header);
   
    // Writes the specified header into the specified BAM file.
    virtual SamStatus::Status writeHeader(IFILE filePtr, 
                                          SamFileHeader& header);

    // Reads the next record from the specified BAM file and stores it in
    // the passed in record.
    virtual void readRecord(IFILE filePtr, 
                            SamFileHeader& header,
                            SamRecord& record, 
                            SamStatus& samStatus);
   
    // Writes the specified record into the specified BAM file.
    virtual SamStatus::Status writeRecord(IFILE filePtr, 
                                          SamFileHeader& header,
                                          SamRecord& record,
                                          SamRecord::SequenceTranslation translation);
};

#endif
