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

#ifndef __SAM_FILE_HEADER_H__
#define __SAM_FILE_HEADER_H__

#include <map>
#include <stdint.h>

#include "SamReferenceInfo.h"
#include "SamHeaderHD.h"
#include "SamHeaderSQ.h"
#include "SamHeaderRG.h"
#include "SamHeaderPG.h"
#include "SamStatus.h"

class SamFileHeader
{
public:
    SamFileHeader();
    ~SamFileHeader();

    // Copy Constructor   
    SamFileHeader(const SamFileHeader& header);

    // Overload operator = to copy the passed in header into this header.
    SamFileHeader & operator = (const SamFileHeader& header);

    // Overload operator = to copy the passed in header into this header.
    bool copy(const SamFileHeader& header);

    void resetHeader();

    // Set the passed in string to the entire header string.  Clearing its
    // current contents.
    // Return true if successfully set (even if set to "")
    bool getHeaderString(std::string& header) const;

    int   getReferenceID(const String & referenceName);
    int   getReferenceID(const char* referenceName);
    const String & getReferenceLabel(int id) const;

    // Get the Reference Information
    const SamReferenceInfo* getReferenceInfo() const;

    // Add reference sequence name and reference sequence length to the header.
    void addReferenceInfo(const char* referenceSequenceName, 
                          int32_t referenceSequenceLength);

    ////////////////////////////////////////////////////////////////////////
    // Set Values in the header
    ////////////////////////////////////////////////////////////////////////

    // Add a header line that is just one tag with a const char* value.
    bool addHeaderLine(const char* type, const char* tag, const char* value); 
    // Add a header line that is already preformatted in a const char*.
    // It is assumed that the line does not contain a \n.
    bool addHeaderLine(const char* headerLine);

//     // Set the specified header type tag to the specified value in the 
//     // header with the specified keyID.  keyID must be specified when
//     // type = SQ, RG, or PG.
//     bool setTag(SamHeaderRecord::SamHeaderRecordType type, const char* tag,
//                 const char* value, const char* keyID = NULL);

    // Set the specified tag to the specified value in the HD header.
    bool setHDTag(const char* tag, const char* value);

    // Set the specified tag to the specified value in the SQ header with
    // the specified name.
    // If the header does not yet exist, the header is added.
    bool setSQTag(const char* tag, const char* value, const char* name);

    // Set the specified tag to the specified value in the RG header with
    // the read group identifier.
    // If the header does not yet exist, the header is added.
    bool setRGTag(const char* tag, const char* value, const char* id);

    // Set the specified tag to the specified value in the PG header with
    // the specified id.
    // If the header does not yet exist, the header is added.
    bool setPGTag(const char* tag, const char* value, const char* id);

    // Add the HD record to the header.
    // Note: it adds a pointer to the passed in header record.  The header
    // record will be deleted when it is cleaned up from this header.
    bool addHD(SamHeaderHD* hd);

    // Add the SQ record to the header.
    // Note: it adds a pointer to the passed in header record.  The header
    // record will be deleted when it is cleaned up from this header.
    bool addSQ(SamHeaderSQ* sq);

    // Add the RG record to the header.
    // Note: it adds a pointer to the passed in header record.  The header
    // record will be deleted when it is cleaned up from this header.
    bool addRG(SamHeaderRG* rg);

    // Add the PG record to the header.
    // Note: it adds a pointer to the passed in header record.  The header
    // record will be deleted when it is cleaned up from this header.
    bool addPG(SamHeaderPG* pg);

    ////////////////////////////////////////////////////////////////////////
    // Remove entries from the header
    ////////////////////////////////////////////////////////////////////////
    bool removeHD();  // Remove the HD record.
    bool removeSQ(const char* name); // Remove SQ record with the specified key.
    bool removeRG(const char* id); // Remove RG record with the specified key.
    bool removePG(const char* id); // Remove PG record with the specified key.


    ////////////////////////////////////////////////////////////////////////
    //
    ////////////////////////////////////////////////////////////////////////
    SamStatus::Status setHeaderFromBamFile(IFILE filePtr);
    
    const char* getHDTagValue(const char* tag);
    // Get the value associated with the specified tag on the SQ line with
    // the specified sequence name.
    const char* getSQTagValue(const char* tag, const char* name);
    // Get the value associated with the specified tag on the RG line with
    // the specified read group identifier.
    const char* getRGTagValue(const char* tag, const char* id);
    // Get the value associated with the specified tag on the RG line with
    // the specified id.
    const char* getPGTagValue(const char* tag, const char* id);

    // Get the number of SQ objects.
    int getNumSQs();

    // Get the number of RG objects.
    int getNumRGs();

    // Get the number of PG objects.
    int getNumPGs();

    // Get the HD object.
    SamHeaderHD* getHD();

    // Get the SQ object with the specified sequence name.
    SamHeaderSQ* getSQ(const char* name);

    // Get the RG object with the specified read group identifier.
    SamHeaderRG* getRG(const char* id);

    // Get the PG object with the specified id.
    SamHeaderPG* getPG(const char* id);

//     //////////////////////////////////
//     // Set methods for header fields.
//     bool setVersion(const char* version);
//     bool setSortOrder(const char* sortOrder);
//     bool addSequenceName(const char* sequenceName);
//     bool setSequenceLength(const char* keyID, int sequenceLength);
//     bool setGenomeAssemblyId(const char* keyID, const char* genomeAssemblyId);
//     bool setMD5Checksum(const char* keyID, const char* md5sum);
//     bool setURI(const char* keyID, const char* uri);
//     bool setSpecies(const char* keyID, const char* species);
//     bool addReadGroupID(const char* readGroupID);
//     bool setSample(const char* keyID, const char* sample);
//     bool setLibrary(const char* keyID, const char* library);
//     bool setDescription(const char* keyID, const char* description);
//     bool setPlatformUnit(const char* keyID, const char* platform);
//     bool setPredictedMedianInsertSize(const char* keyID, const char* isize);
//     bool setSequencingCenter(const char* keyID, const char* center);
//     bool setRunDate(const char* keyID, const char* runDate);
//     bool setTechnology(const char* keyID, const char* technology);
//     bool addProgram(const char* programID);
//     bool setProgramVersion(const char* keyID, const char* version);
//     bool setCommandLine(const char* keyID, const char* commandLine);
    
//     ///////////////////////////////////
//     // Get methods for header fields.
//     // Returns the number of SQ entries in the header.
//     int32_t getSequenceDictionaryCount();
    // Return the Sort Order value that is set in the Header.
    // If this field does not exist, "" is returned.
    const char* getSortOrder();


    // DEPRECATED
    const char* getTagSO();

    // Get the next SQ header record.  After all SQ headers have been retrieved,
    // NULL is returned until a reset is called.
    SamHeaderRecord* getNextSQRecord();

    // Get the next RG header record.  After all RG headers have been retrieved,
    // NULL is returned until a reset is called.
    SamHeaderRecord* getNextRGRecord();

    // Get the next PG header record.  After all PG headers have been retrieved,
    // NULL is returned until a reset is called.
    SamHeaderRecord* getNextPGRecord();

    // Reset to the beginning of the header records so the next call
    // to getNextSQRecord returns the first SQ header record.
    void resetSQRecordIter();

    // Reset to the beginning of the header records so the next call
    // to getNextRGRecord returns the first RG header record.
    void resetRGRecordIter();

    // Reset to the beginning of the header records so the next call
    // to getNextPGRecord returns the first PG header record.
    void resetPGRecordIter();

    // Get the next header record of the specified type.
    // Pass in the index to start looking at and the type to look for.
    // Update the index.
    // After all headers of that type have been retrieved,
    // NULL is returned until a reset is called for that type.
    SamHeaderRecord* getNextHeaderRecord(uint32_t& index, 
                                         SamHeaderRecord::SamHeaderRecordType headerType);

    // Get the next header record.  After all headers have been retrieved,
    // NULL is returned until a reset is called.  Does not return the
    // Comment lines.
    // NOTE: both getNextHeaderRecord and getNextHeaderLine increment the
    // same iterator.
    SamHeaderRecord* getNextHeaderRecord();


    // Set the passed in string to the next header line.  The passed in 
    // string will be overwritten.  If there are no more header lines or there
    // is an error, false is returned and the passed in string is set to ""
    // until a rest is called.
    // Will also return the comment lines.
    // NOTE: both getNextHeaderRecord and getNextHeaderLine increment the
    // same iterator.
    bool getNextHeaderLine(std::string &headerLine);

    // Reset to the beginning of the header records so the next call
    // to getNextHeaderRecord returns the first header line.
    void resetHeaderRecordIter();
   
    // Returns the comment on the next comment line.  Returns "" if all comment
    // lines have been returned, until resetCommentIter is called.
    const char* getNextComment();

    // Resets to the beginning of the comments so getNextComment returns
    // the first comment.
    void resetCommentIter();

    // Add a comment.
    bool addComment(const char* comment);

    // Populate the reference info from the SQ fields.
    void generateReferenceInfo();


private:
    // Parse the header string.
    bool parseHeader(String& header);

    // Parse the specified line of the header.
    bool parseHeaderLine(const String& headerLine);

    // Set the passed in string to the header line at the specified index.
    // It does NOT clear the current contents of header.
    bool getHeaderLine(unsigned int index, std::string& header) const;

    int16_t makeKey(char ch1, char ch2)
    {
        return((ch1 << 8) + ch2);
    }

    // Only one HD type is allowed per file.
    SamHeaderHD* myHD;

    // There can be multiple SQ Types, indexed by SN.
    StringHash mySQs;

    // There can be multiple RG Types, indexed by ID.
    StringHash myRGs;

    // There can be multiple PG types, indexed by ID.
    StringHash myPGs;

    // Reference Name information
    SamReferenceInfo myReferenceInfo;

    // Vector of comments
    std::vector<std::string> myComments;

    std::vector<SamHeaderRecord*> myHeaderRecords;

    uint32_t myCurrentSQIndex;

    uint32_t myCurrentRGIndex;

    uint32_t myCurrentPGIndex;

    uint32_t myCurrentHeaderIndex;

    uint32_t myCurrentCommentIndex;

    static const std::string EMPTY_RETURN;
};

#endif

