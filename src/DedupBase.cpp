/*
 *  Copyright (C) 2010-2016  Regents of the University of Michigan
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
#include "DedupBase.h"
#include "Logger.h"
#include "SamFlag.h"

DedupBase::~DedupBase()
{
    mySecondarySupplementaryMap.clear();
}

void DedupBase::markDuplicateLoop(bool verboseFlag, bool removeFlag, const String& inFile, const String& outFile, Recab* recabPtr)
{
    static bool firstMapUnmapOrderError = true;
    ////////////////////////////////////////////////////////////
    // Second Pass through the file
    SamFile samIn;
    SamFileHeader header;
    SamFile samOut;

    // get ready to write the output file by making a second pass
    // through the input file
    samIn.OpenForRead(inFile.c_str());
    samIn.ReadHeader(header);

    samOut.OpenForWrite(outFile.c_str());
    samOut.WriteHeader(header);

    // If we are recalibrating, output the model information.
    if(recabPtr != NULL)
    {
        recabPtr->modelFitPrediction(outFile);
    }
    // an iterator to run through the duplicate indices
    int currentDupIndex = 0;
    bool moreDups = !myDupList.empty();

    // let the user know what we're doing
    Logger::gLogger->writeLog("\nWriting %s", outFile.c_str());
    
    // count the duplicate records as a check
    uint32_t singleDuplicates(0), pairedDuplicates(0),
        mateUnmappedDuplicates(0),
        secondaryDuplicates(0), supplementaryDuplicates(0);

    // start reading records and writing them out
    SamRecord *prevRecord = NULL;
    SamRecord *recordPtr = mySamPool.getRecord();
    while(samIn.ReadRecord(header, *recordPtr))
    {
        uint32_t currentIndex = samIn.GetCurrentRecordCount();

        bool foundDup = moreDups &&
            (currentIndex == myDupList[currentDupIndex]);

        // modify the duplicate flag and write out the record,
        // if it's appropriate
        int flag = recordPtr->getFlag();

        if (foundDup)
        {   
            currentDupIndex++;
            // increment duplicate counters to verify we found them all
            if ( ( flag & 0x0001 ) == 0 )
            { // unpaired
                singleDuplicates++;
            }
            else if ( flag & 0x0008 )
            { // mate unmapped
                mateUnmappedDuplicates++;
            }
            else
            {
                pairedDuplicates++;
            }

            // Update NonPrimaryMap to indicate this readname is a duplicate
            // (if the readname is in the map)
            // Need to do it here in case the non-primary read comes after the primary one
            // so the read name was not in the map when the primary read was encountered
            // on the first pass through the file.
            if(myMarkNonPrimary)
            {
                markDuplicateInNonPrimaryMaps(recordPtr);
            }
        }
        else if(myMarkNonPrimary && (SamFlag::isSecondary(flag) || 
                                     (flag & SamFlag::SUPPLEMENTARY_ALIGNMENT)))
        {
            // Check the supplemenatry map to see if this read is a duplicate.
            foundDup = mySecondarySupplementaryMap[recordPtr->getReadName()];
            if(foundDup)
            {
                if(flag & SamFlag::SUPPLEMENTARY_ALIGNMENT)
                {
                    ++supplementaryDuplicates;
                }
                else
                {
                    ++secondaryDuplicates;
                }
            }
        }

        setDuplicate(recordPtr, foundDup);
        // Update the flag value
        flag = recordPtr->getFlag();

        // Handle previous record if there is one
        if(prevRecord != NULL)
        {
            if(strcmp(prevRecord->getReadName(), recordPtr->getReadName()) == 0)
            {
                // Same Read Name, so update the two mates to have the same duplicate marking as the
                // mapped mate.
                if(SamFlag::isMapped(flag))
                {
                    // Current record is mapped, so prevRecord shouldn't be - copy the flag from
                    // the current record to the previous one.
                    setDuplicate(prevRecord, SamFlag::isDuplicate(flag));
                }
                else
                {
                    // Current record is not mapped, so prevRecord should be - copy the flag from
                    // the previous record to the current one.
                    setDuplicate(recordPtr, SamFlag::isDuplicate(prevRecord->getFlag()));
                    // Update the flag value
                    flag = recordPtr->getFlag();
                }
                // recalibrate prevRecord if necessary.
                if(recabPtr != NULL)
                {
                    recabPtr->processReadApplyTable(*prevRecord);
                }
                // Write the previousRecord
                if(!SamFlag::isDuplicate(flag) || (!removeFlag ))
                {
                    samOut.WriteRecord(header, *prevRecord);
                }
                // Clear prevRecord since it has been handled.
                mySamPool.releaseRecord(prevRecord);
                prevRecord = NULL;
                // Done processing the pair
            }
            else
            {
                // There is a previous record, but it has a different read name, so write a warning
                if(firstMapUnmapOrderError)
                {
                    Logger::gLogger->warning("Unmapped mate not found next to the mapped mate for readname %s. Additional warnings for this are suppressed.",
                                             prevRecord->getReadName());
                    std::cerr << "WARNING: Unmapped mate not found next to the mapped mate for readname " <<
                        prevRecord->getReadName() << ". Additional warnings for this are suppressed.\n";
                    firstMapUnmapOrderError = false;
                }
                // recalibrate prevRecord if necessary.
                if(recabPtr != NULL)
                {
                    recabPtr->processReadApplyTable(*prevRecord);
                }
                // Write the previous record if it is not a duplicate or  if we are not removing duplicates
                if(!SamFlag::isDuplicate(prevRecord->getFlag()) || (!removeFlag ))
                {
                    samOut.WriteRecord(header, *prevRecord);
                }
                // Clear prevRecord since it has been handled.
                mySamPool.releaseRecord(prevRecord);
                prevRecord = NULL;

                // Check if the current record needs to be stored.
                if(SamFlag::isMapped(flag) != SamFlag::isMateMapped(flag))
                {
                    // Current record & mate - one is mapped, one is not.
                    // No previous record, so the mate must be next.  Store this record so the next read can
                    // be checked for the mate.
                    prevRecord = recordPtr;
                }
            }
        }
        else
        {      
            // No previous record.  Check if the current record needs to be stored.
            if(SamFlag::isMapped(flag) != SamFlag::isMateMapped(flag))
            {
                // Current record & mate - one is mapped, one is not.
                // No previous record, so the mate must be next.  Store this record so the next read can
                // be checked for the mate.
                prevRecord = recordPtr;
            }
        }
        
        // write the record if we didn't save it as the previous record, and
        // we are not removing duplicates
        if(prevRecord == recordPtr)
        {
            // Stored this record, so get a new one for the next read.
            recordPtr = mySamPool.getRecord();
        }
        else
        {
            // Didn't store the read for processing with the next read
            // recalibrate if necessary.
            if(recabPtr != NULL)
            {
                recabPtr->processReadApplyTable(*recordPtr);
            }
            if (!foundDup || !removeFlag ) samOut.WriteRecord(header, *recordPtr);
        }

        // Let the user know we're still here
        if (verboseFlag && (currentIndex % 100000 == 0)) {
            Logger::gLogger->writeLog("recordCount=%u", currentIndex);
        }
    }

    // We're done.  Close the files and print triumphant messages.
    samIn.Close();
    samOut.Close();

    Logger::gLogger->writeLog("Successfully %s %u unpaired, %u mate unmapped, %u paired duplicate reads, %u secondary reads, and %u supplementary reads", 
                              removeFlag ? "removed" : "marked" ,
                              singleDuplicates,
                              mateUnmappedDuplicates,
                              pairedDuplicates/2,
                              secondaryDuplicates,
                              supplementaryDuplicates);
}


void DedupBase::setDuplicate(SamRecord* recordPtr, bool duplicate)
{
    if(duplicate)
    {
        recordPtr->setFlag( recordPtr->getFlag() | 0x400 );
    }
    else if(myForceFlag)
    {
        // this is not a duplicate we've identified but we want to
        // remove any duplicate marking
        recordPtr->setFlag( recordPtr->getFlag() & 0xfffffbff ); // unmark duplicate
    }
}


void DedupBase::markDuplicateInNonPrimaryMaps(SamRecord* recordPtr)
{
    // Check the NonPrimary Maps.
    if(recordPtr == NULL)
    {
        return;
    }
    int flag = recordPtr->getFlag();

    // check the maps if the read is:
    //     * single end
    //     * the mate is unmapped (treated as single end)
    //  or * first fragment (only want to check for 1 of the 2 mates)
    if(!SamFlag::isPaired(flag) || !SamFlag::isMateMapped(flag) ||
       SamFlag::isFirstFragment(flag))
    {
        NonPrimaryMap::iterator iter;
        iter = mySecondarySupplementaryMap.find(recordPtr->getReadName());
        if(iter != mySecondarySupplementaryMap.end())
        {
            // Found in the non primary map, so set it to be a dupilcate.
            iter->second = true;
        }
    }
}
