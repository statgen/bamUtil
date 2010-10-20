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

#include "SamStatus.h"

const char* SamStatus::enumStatusString[] = {
    "SUCCESS",
    "UNKNOWN",
    "NO_MORE_RECS",
    "FAIL_IO",
    "FAIL_ORDER",
    "FAIL_PARSE",
    "INVALID_SORT",
    "INVALID", 
    "FAIL_MEM"
};


const char* SamStatus::getStatusString(SamStatus::Status statusEnum)
{
    return(enumStatusString[statusEnum]);
}


// Returns whether or not it is "safe" to keep processing the file
// after the specified status return.
bool SamStatus::isContinuableStatus(SamStatus::Status status)
{
    if(status == SamStatus::SUCCESS || status == SamStatus::FAIL_PARSE || 
       status == SamStatus::INVALID_SORT || status == SamStatus::INVALID)
    {
        // The status is such that file processing can continue.
        return(true);
    }
    // UNKNOWN, NO_MORE_RECS, FAIL_IO, FAIL_ORDER, FAIL_MEM
    return(false);
}


// Constructor
SamStatus::SamStatus(ErrorHandler::HandlingType handleType)
    : myHandlingType(handleType)
{
    reset();
}

   
// Destructor
SamStatus::~SamStatus()
{
}


// Resets this status.
void SamStatus::reset()
{
    myType = UNKNOWN;
    myMessage = "";
}


void SamStatus::setHandlingType(ErrorHandler::HandlingType handleType)
{
    myHandlingType = handleType;
}


// Set the status with the specified values.
void SamStatus::setStatus(Status newStatus, const char* newMessage)
{
    myType = newStatus;
    myMessage = getStatusString(newStatus);
    myMessage += ": ";
    myMessage += newMessage;

    if(newStatus != SUCCESS)
    {
        handleError(newStatus, newMessage);
    }
}


// Adds the specified error message to the status message.
// Sets the status to newStatus if the current status is SUCCESS.
void SamStatus::addError(Status newStatus, const char* newMessage)
{
    if(myType == SamStatus::SUCCESS)
    {
        myType = newStatus;
    }
    else
    {
        myMessage += "\n";
    }
    myMessage += getStatusString(newStatus);
    myMessage += ": ";
    myMessage += newMessage;

    if(newStatus != SUCCESS)
    {
        handleError(newStatus, newMessage);
    }
}


// Adds the specified status to the status message.
// Sets the status to newStatus if the current status is SUCCESS.
void SamStatus::addError(SamStatus newStatus)
{
    if(myType == SamStatus::SUCCESS)
    {
        myType = newStatus.myType;
    }
    else
    {
        myMessage += "\n";
    }
    myMessage += newStatus.myMessage;

    if(newStatus != SUCCESS)
    {
        handleError(newStatus.myType, newStatus.myMessage.c_str());
    }
}


// Return the enum for this status.
SamStatus::Status SamStatus::getStatus() const
{
    return(myType);
}


// Return the status message.
const char* SamStatus::getStatusMessage() const
{
    return(myMessage.c_str());
}


// Overload operator = to set the sam status type to the
// passed in status and to clear the message string.
SamStatus & SamStatus::operator = (SamStatus::Status newStatus)
{
    reset();
    myType = newStatus;

    if(newStatus != SUCCESS)
    {
        handleError(newStatus, "");
    }
    return(*this);
}


// Overload operator != to determine if the passed in type is not equal
// to this status's type.
bool SamStatus::operator != (const SamStatus::Status& compStatus) const
{
    return(compStatus != myType);
}


// Overload operator != to determine if the passed in type is equal
// to this status's type.
bool SamStatus::operator == (const SamStatus::Status& compStatus) const
{
    return(compStatus == myType);
}


void SamStatus::handleError(Status newStatus, const char* newMessage)
{
    // If the status is not success and not NO_MORE_RECS, handle
    // the error  (SUCCESS & NO_MORE_RECS are not real errors.)
    if((newStatus != SUCCESS) && (newStatus != NO_MORE_RECS))
    {
        std::string message = getStatusString(newStatus);
        message += ": ";
        message += newMessage;
 
        ErrorHandler::handleError(message.c_str(), myHandlingType);
    }
}
