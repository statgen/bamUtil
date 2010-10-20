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

#ifndef __SAM_VALIDATION_H__
#define __SAM_VALIDATION_H__

#include "SamFile.h"
#include <list>

class SamValidationError
{
public:
    // Warning is used if it is just an invalid value.
    // Error is used if parsing could not succeed.
    enum Severity {WARNING, ERROR};

    enum Type
        {
            // 
            INVALID_QNAME,
            INVALID_REF_ID,
            INVALID_RNAME,
            INVALID_POS,
            INVALID_MAPQ,
            INVALID_CIGAR,
            INVALID_MRNM,
            INVALID_QUAL
        };

    static const char* getTypeString(Type type);

    SamValidationError(Type type, Severity severity, std::string Message);
   
    Type getType() const;
    Severity getSeverity() const;
    const char* getMessage() const;

    const char* getTypeString() const;
    const char* getSeverityString() const;

    void getErrorString(std::string& errorString) const;

    void printError() const;

private:
    SamValidationError();

    static const char* enumTypeString[];
    static const char* enumSeverityString[];

    Type myType;
    Severity mySeverity;
    std::string myMessage;

};


//
// stream output for validation failure information
//
inline std::ostream &operator << (std::ostream &stream, 
                                  const SamValidationError &error)
{
    std::string errorMessage;
    error.getErrorString(errorMessage);
    stream << errorMessage;
    return stream;
}


class SamValidationErrors
{
public:
    // Constructor.
    SamValidationErrors();
    // Destructor
    ~SamValidationErrors();

    // Remove all the errors from the list.
    void clear();

    // Adds the specified error to the set of errors.
    void addError(SamValidationError::Type newType, 
                  SamValidationError::Severity newSeverity,
                  const char* newMessage);

    // Return the number of validation errors that are contained in this object.
    unsigned int numErrors();

    // Return a pointer to the next error.  It does not remove it from the list.
    // Returns null once all errors have been retrieved until resetErrorIter
    // is called.
    const SamValidationError* getNextError();
   
    // Resets the iterator to the begining of the errors.
    void resetErrorIter();

    // Appends the error messages to the passed in string.
    void getErrorString(std::string& errorString) const;

private:
    std::list<const SamValidationError*> myValidationErrors;
    std::list<const SamValidationError*>::const_iterator myErrorIter;
};


//
// stream output for all validation failures information
//
inline std::ostream& operator << (std::ostream& stream,
                                  const SamValidationErrors& errors)
{
    std::string errorString = "";
    errors.getErrorString(errorString);
    stream << errorString;
    return stream;
}


class SamValidator
{
public:

    static bool isValid(SamFileHeader& samHeader, SamRecord& samRecord, 
                        SamValidationErrors& validationErrors);

    static bool isValidQname(const char* qname, uint8_t qnameLen, 
                             SamValidationErrors& validationErrors);
    static bool isValidFlag(uint16_t flag,
                            SamValidationErrors& validationErrors);
    // Validate the rname including validating against the header.
    static bool isValidRname(SamFileHeader& samHeader, 
                             const char* rname,
                             SamValidationErrors& validationErrors);
    // Validate the rname without validating against the header.
    static bool isValidRname(const char* rname,
                             SamValidationErrors& validationErrors);
    static bool isValidRefID(int32_t refID, const StringArray& refContigs, 
                             SamValidationErrors& validationErrors);
    static bool isValid1BasedPos(int32_t pos, 
                                 SamValidationErrors& validationErrors);
    static bool isValidMapQuality(uint8_t mapQuality,
                                  SamValidationErrors& validationErrors);
    // Cigar validation depends on sequence.
    static bool isValidCigar(const char* cigar, const char* sequence,
                             SamValidationErrors& validationErrors);
    static bool isValidMrnm();
    static bool isValidMpos();
    static bool isValidIsize();
    static bool isValidSeq();
    // Quality validation depends on sequence.
    static bool isValidQuality(const char* quality, const char* sequence,
                                  SamValidationErrors& validationErrors);
    static bool isValidTag();
    static bool isValidVtype();
    static bool isValidValue();
};


#endif
