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

#ifndef __SAMHEADER_RECORD_H__
#define __SAMHEADER_RECORD_H__

#include "StringArray.h"
#include "StringHash.h"
#include "SamHeaderTag.h"

class SamHeaderRecord
{
public:
    enum SamHeaderRecordType {HD, SQ, RG, PG};

    // Constructor
    SamHeaderRecord();
   
    // Destructor
    ~SamHeaderRecord();

    // Set the fields from the passed in line.
    // Return true if successfully set.
    bool setFields(const StringArray& tokens);

    // Check to see if the record is valid.
    bool isValid();

    // Return the value associated with the specified tag.
    const char* getTagValue(const char* tag) const;

    // Set the value of the specified tag to the specified value.
    // Set value to NULL in order to delete the tag.
    // Returns whether or not it was successful.
    // Fails if tag is the key tag and the key tag already exists.
    bool setTag(const char* tag, const char* value);

    // Reset this header record to an empty state.
    void reset();

    // Appends the string representation of this header record
    // to the passed in string.
    bool appendString(std::string& header);

    // Add the key tag with the specified value.
    bool addKey(const char* value);

    // This record is active if there is at least one tag set.
    bool isActiveHeaderRecord();

    // Return the type of this header record.
    const char* getTypeString();

    // Return the type of this header record.
    SamHeaderRecordType getType();

protected:
    void addRequiredTag(const char* requiredTag);

    // The type for this header record.
    std::string myTypeString;

    // The type for this header record.
    SamHeaderRecordType myType;

    // The TAG name that is the key for this record
    // Only applicable if more than one of this type
    // of record is allowed.
    std::string myKeyTag;

private:
    // hash from tag name to index into the tag values vector.
    StringIntHash myTagHash;
    std::vector<SamHeaderTag*> myTags;

    // The tags that are required for this record.
    std::vector<String> myRequiredTags;

    int myNumActiveTags;
};

#endif
