/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

#include "BaseInfoFile.h"

BaseInfoFileReader::BaseInfoFileReader()
{
}


BaseInfoFileReader::~BaseInfoFileReader()
{
    close();
}


bool BaseInfoFileReader::open(const char* fileName)
{
    // Close any already open file.
    close();
    
    myFilePtr = ifopen(filename, "r");
    if(myFilePtr = NULL)
    {
        std::string errorMessage = "Failed to Open ";
        errorMessage += filename;
        errorMessage += " for reading.";
        myStatus.setStatus(StatGenStatus::FAIL_IO, errorMessage.c_str());
        return(false);
    }
    return(true);
}


void BaseInfoFileReader::close()
{
    ifclose(myFilePtr);
}


bool BaseInfoFileReader::isOpen()
{
    if (myFilePtr != NULL)
    {
        // File Pointer is set, so return if it is open.
        return(myFilePtr->isOpen());
    }
    // File pointer is not set, so return false, not open.
    return false;
}


bool BaseInfoFileReader::isEof()
{
    if (myFilePtr != NULL)
    {
        // File Pointer is set, so return if eof.
        return(ifeof(myFilePtr));
    }
    // File pointer is not set, so return true, eof.
    return true;
}


bool BaseInfoFileReader::getNextDataRecord(BaseInfoRecord& rec)
{
    
}
