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

#ifndef __SAM_STATUS_H__
#define __SAM_STATUS_H__

#include <iostream>
#include "ErrorHandler.h"

class SamStatus
{
public:

    // Return value enum for the SamFile class methods.
    //    SUCCESS      : method completed successfully.
    //    UNKNOWN      : unknown result (default value should never be used)
    //    NO_MORE_RECS : failed to read a record since there are no more to read
    //                   either in the file or section if section based reading.
    //    FAIL_IO      : method failed due to an I/O issue.
    //    FAIL_ORDER   : method failed because it was called out of order,
    //                   like trying to read a file without opening it for
    //                   read or trying to read a record before the header.
    //    FAIL_PARSE   : failed to parse a record/header - invalid format.
    //    INVALID_SORT : record is invalid due to it not being sorted.
    //    INVALID      : record is invalid other than for sorting.
    //    FAIL_MEM     : fail a memory allocation.
    enum Status {SUCCESS = 0, UNKNOWN, NO_MORE_RECS, FAIL_IO, FAIL_ORDER,
                 FAIL_PARSE, INVALID_SORT, INVALID, FAIL_MEM};

    static const char* getStatusString(SamStatus::Status statusEnum);

    // Returns whether or not it is "safe" to keep processing the file
    // after the specified status return.
    static bool isContinuableStatus(SamStatus::Status status);

    // Constructor
    SamStatus(ErrorHandler::HandlingType handleType = ErrorHandler::EXCEPTION);
   
    // Destructor
    ~SamStatus();

    // Resets this status.
    void reset();

    void setHandlingType(ErrorHandler::HandlingType handleType);

    // Set the status with the specified values.
    void setStatus(Status newStatus, const char* newMessage);

    // Adds the specified error message to the status message.
    // Sets the status to newStatus if the current status is SUCCESS.
    void addError(Status newStatus, const char* newMessage);


    // Adds the specified status to the status message.
    // Sets the status to newStatus if the current status is SUCCESS.
    void addError(SamStatus newStatus);

    // Return the enum for this status.
    Status getStatus() const;

    // Return the status message.
    const char* getStatusMessage() const;

    // Overload operator = to set the sam status type to the
    // passed in status and to clear the message string.
    SamStatus & operator = (Status newStatus);
   
    //    // Overload operator = to set the sam status.
    //    SamStatus & operator = (SamStatus newStatus);
   
    // Overload operator != to determine if the passed in type is not equal
    // to this status's type.
    bool operator != (const SamStatus::Status& compStatus) const;

    // Overload operator != to determine if the passed in type is equal
    // to this status's type.
    bool operator == (const SamStatus::Status& compStatus) const;
      
private:
    // Handle an error based on the error handling type.
    void handleError(Status newType, const char* newMessage);


    static const char* enumStatusString[];

    Status myType;
    std::string myMessage;
    ErrorHandler::HandlingType myHandlingType;
};


#endif
