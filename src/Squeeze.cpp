/*
 *  Copyright (C) 2010-2015  Regents of the University of Michigan
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

#include "Squeeze.h"
#include "BgzfFileType.h"
// #include <stdio.h>
// #include <string.h>
// #include "SamFile.h"
#include "SamFlag.h"

Squeeze::Squeeze()
    : myBinMid(false),
      myBinHigh(false),
      myBinQualS(""),
      myBinQualF("")
{
    for(int i = 0; i <= MAX_QUAL_CHAR; i++)
    {
        myQualBinMap[i] = 0;
    }
}

Squeeze::~Squeeze()
{
}

void Squeeze::squeezeDescription()
{
    std::cerr << " squeeze -  reduces files size by dropping OQ fields, duplicates, & specified tags, using '=' when a base matches the reference, binning quality scores, and replacing readNames with unique integers" << std::endl;
}


void Squeeze::description()
{
    squeezeDescription();
}


void Squeeze::binningUsageLine()
{
    std::cerr << "[--binQualS <minQualBin2>,<minQualBin3><...>] [--binQualF <filename>] [--binMid|binHigh]";
}


void Squeeze::binningUsage()
{
    std::cerr << "\tQuality Binning Parameters (optional):" << std::endl;
    std::cerr << "\t  Bin qualities by phred score, into the ranges specified by binQualS or binQualF (both cannot be used)" << std::endl;
    std::cerr << "\t  Ranges are specified by comma separated minimum phred score for the bin, example: 1,17,20,30,40,50,70" << std::endl;
    std::cerr << "\t  The first bin always starts at 0, so does not need to be specified." << std::endl;
    std::cerr << "\t  By default, the bin value is the low end of the range." << std::endl;
    std::cerr << "\t\t--binQualS   : Bin the Qualities as specified (phred): minQualOfBin2, minQualofBin3..." << std::endl;
    std::cerr << "\t\t--binQualF   : Bin the Qualities based on the specified file" << std::endl;
    std::cerr << "\t\t--binMid     : Use the mid point of the quality bin range for the quality value of the bin." << std::endl;
    std::cerr << "\t\t--binHigh    : Use the high end of the quality bin range for the quality value of the bin." << std::endl;
}


// print Usage
void Squeeze::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam squeeze --in <inputFile> --out <outputFile.sam/bam/ubam (ubam is uncompressed bam)> [--refFile <refFilePath/Name>] [--keepOQ] [--keepDups] [--readName <readNameMapFile.txt>] [--sReadName <readNameMapFile.txt>] [--rmTags <Tag:Type[,Tag:Type]*>] [--noeof] [--params] ";
    binningUsageLine();
    std::cerr << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in         : the SAM/BAM file to be read" << std::endl;
    std::cerr << "\t\t--out        : the SAM/BAM file to be written" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--refFile    : reference file name used to convert any bases that match the reference to '='" << std::endl;
    std::cerr << "\t\t--keepOQ     : keep the OQ tag rather than removing it.  Default is to remove it." << std::endl;
    std::cerr << "\t\t--keepDups   : keep duplicates rather than removing records marked duplicate.  Default is to remove them." << std::endl;
    std::cerr << "\t\t--sReadName  : Replace read names with unique integers and write the mapping to the specified file." << std::endl;
    std::cerr << "                   This version requires the input file to have been presorted by readname, but" << std::endl;
    std::cerr << "                   no validation is done to ensure this.  If it is not sorted, a readname will" << std::endl;
    std::cerr << "                   get mapped to multiple new values." << std::endl;
    std::cerr << "\t\t--readName   : Replace read names with unique integers and write the mapping to the specified file." << std::endl; 
    std::cerr << "                   This version does not require the input file to have been presorted by readname," << std::endl;
    std::cerr << "                   but uses a lot of memory since it stores all the read names." << std::endl;
    std::cerr << "\t\t--rmTags     : Remove the specified Tags formatted as Tag:Type,Tag:Type,Tag:Type..." << std::endl;
    std::cerr << "\t\t--noeof      : do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params     : print the parameter settings" << std::endl;
    binningUsage();
    std::cerr << std::endl;
}

// main function
int Squeeze::execute(int argc, char ** argv)
{
    String inFile = "";
    String outFile = "";
    String refFile = "";
    bool noeof = false;
    bool params = false;
    bool keepOQ = false;
    bool keepDups =  false;
    String readName = "";
    String sReadName = "";
    IFILE readNameFile = NULL;
    myBinMid = false;
    myBinHigh = false;
    String rmTags = "";

    ParameterList inputParameters;
    LongParamContainer parameters;

    parameters.addGroup("Required Parameters");
    parameters.addString("in", &inFile);
    parameters.addString("out", &outFile);
    parameters.addGroup("Optional Parameters");
    parameters.addString("refFile", &refFile);
    parameters.addBool("keepOQ", &keepOQ);
    parameters.addBool("keepDups", &keepDups);
    parameters.addString("readName", &readName);
    parameters.addString("sReadName", &sReadName);
    parameters.addString("rmTags", &rmTags);
    parameters.addBool("noeof", &noeof);
    parameters.addBool("params", &params);
    parameters.addPhoneHome(VERSION);
    addBinningParameters(parameters);    

    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            parameters.getLongParameterList()));
    
    // parameters start at index 2 rather than 1.
    inputParameters.Read(argc, argv, 2);
    
    // If no eof block is required for a bgzf file, set the bgzf file type to 
    // not look for it.
    if(noeof)
    {
        // Set that the eof block is not required.
        BgzfFileType::setRequireEofBlock(false);
    }
    
    // Check to see if the in file was specified, if not, report an error.
    if(inFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(outFile == "")
    {
        usage();
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Check that both readName and sortedRN aren't specified since
    // they mean the same thing, except the one indicates that the input file is sorted by read name.
    if(!readName.IsEmpty() && !sReadName.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "ERROR: --readName and --sReadName cannot both be specified\n";
        return(-1);
    }

    int status = processBinningParam();
    if(status != 0)
    {
        std::cerr << "HI\n";
        inputParameters.Status();
        return(status);
    }

    // Setup the read name map file.
    if(!sReadName.IsEmpty())
    {
        readName = sReadName;
    }
    if(!readName.IsEmpty())
    {
        readNameFile = ifopen(readName, "w");
        if(readNameFile == NULL)
        {
            std::cerr << "Failed to open the readName File for write: " << readName << std::endl;
            return(-1);
        }
    }

    if(params)
    {
        inputParameters.Status();
    }

    // Open the input file for reading.
    SamFile samIn;
    samIn.OpenForRead(inFile);

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outFile);
    // Check to see if the ref file was specified.
    // Open the reference.
    GenomeSequence* refPtr = NULL;
    if(refFile != "")
    {
        refPtr = new GenomeSequence(refFile);
        // Since a reference was specified, convert matching bases to '='.
        samOut.SetWriteSequenceTranslation(SamRecord::EQUAL);
        // Set the reference for the output file.
        samOut.SetReference(refPtr);
    }

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    // Write the sam header.
    samOut.WriteHeader(samHeader);

    SamRecord samRecord;

    // Create the hash for readnames.
    StringIntHash rnHash;
    // Track the next read name to assign.
    int32_t nextRn = 0;
    int32_t numHashEntries = rnHash.Entries();
    String newRn = "";
    String prevRn = "";

    // Set returnStatus to success.  It will be changed to the
    // failure reason if any of the writes or updates fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;
  
    // Keep reading records until ReadRecord returns false.
    while(samIn.ReadRecord(samHeader, samRecord))
    {
        // Successfully read a record from the file.

        // Remove the record if it is a duplicate and we are not
        // supposed to keep duplicates.
        if(!keepDups && SamFlag::isDuplicate(samRecord.getFlag()))
        {
            // Duplicate, so do not write it to the output
            // file and just continue to the next record.
            continue;
        }

        if(readNameFile != NULL)
        {
            // Shorten the readname.
            // Check the hash for the readname.
            const char* readName = samRecord.getReadName();
            
            if(sReadName.IsEmpty())
            {
                // Lookup the readname in the hash.
                int index = rnHash.Find(readName, nextRn);
                
                // Check to see if the nextRn was added to the hash.
                if(rnHash.Entries() != numHashEntries)
                {
                    // New Read Name was added to the hash.
                    numHashEntries = rnHash.Entries();
                    newRn = nextRn;
                    
                    // Write it to the file.
                    ifprintf(readNameFile, "%s\t%d\n", readName, nextRn);
                    
                    // Update the next read name.
                    ++nextRn;
                }
                else
                {
                    // Found the read name, so use that value.
                    newRn = rnHash.Integer(index);
                }
            }
            else
            {
                // the read is sorted, so check if it is the same as the previous
                // read name.
                if(prevRn != readName)
                {
                    // New read name 
                    newRn = nextRn;
                    
                    // Write it to the file.
                    ifprintf(readNameFile, "%s\t%d\n", readName, nextRn);
                    
                    // Update the next read name.
                    ++nextRn;
                    prevRn = readName;
                }
                
            }
            samRecord.setReadName(newRn.c_str());
        }

        // Remove the OQ tag if we are not supposed to remove OQ tags.
        if(!keepOQ)
        {
            if(!samRecord.rmTag("OQ", 'Z'))
            {
                // Failed to remove a tag.
                SamStatus errorStatus = samRecord.getStatus();
                fprintf(stderr, "%s\n", errorStatus.getStatusMessage());
                returnStatus = errorStatus.getStatus();
            }
        }

        // Remove any specified tags.
        if(!rmTags.IsEmpty())
        {
            if(!samRecord.rmTags(rmTags.c_str()))
            {
                // Failed to remove the specified tags.
                fprintf(stderr, "%s\n", samIn.GetStatusMessage());
                returnStatus = samIn.GetStatus();
            }
        }

        // Bin the qualities.
        bin(samRecord);

        if(!samOut.WriteRecord(samHeader, samRecord))
        {
            // Failed to write a record.
            fprintf(stderr, "Failure in writing record %s\n", samOut.GetStatusMessage());
            returnStatus = samOut.GetStatus();
        }
    }
   
    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        returnStatus = samOut.GetStatus();
    }   
   
    std::cerr << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cerr << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    // Since the reads were successful, return the status based
    // on the status of the reads/writes.  If any failed, return
    // their failure status.
    samIn.Close();
    samOut.Close();
    return returnStatus;
}


void Squeeze::addBinningParameters(LongParamContainer& params)
{
    params.addGroup("Optional Quality Binning Parameters");
    params.addString("binQualS", &myBinQualS);
    params.addString("binQualF", &myBinQualF);
    params.addBool("binMid", &myBinMid);
    params.addBool("binHigh", &myBinHigh);
}


int Squeeze::processBinningParam()
{
    if(!myBinQualS.IsEmpty() && !myBinQualF.IsEmpty())
    {
        std::cerr << "ERROR: --binQualS and --binQualF cannot both be specified\n";
        return(-1);
    }

    if(!myBinQualF.IsEmpty())
    {
        // Read the quality bins from the file.
        IFILE qualBinFile = ifopen(myBinQualF, "r");
        if(qualBinFile == NULL)
        {
            std::cerr << "ERROR: failed to open the quality bin file (" 
                      << qualBinFile << ")." << std::endl;
            return(-1);
        }

        if(myBinQualS.ReadLine(qualBinFile) <= 0)
        {
            std::cerr << "ERROR: failed to read the quality bin file (" 
                      << qualBinFile << ")." << std::endl;
            return(-1);
        }
    }

    if(!myBinQualS.IsEmpty())
    {
        // Determine the bins.
        StringArray bins;
        bins.ReplaceColumns(myBinQualS, ',');
        int binStart = 0;
        int nextBinStart = 0;

        // Fill the bins, by reading the starting bin positions.
        // The previous bin ends just before the next one starts.
        for(int i = 0; i < bins.Length(); i++)
        {
            // The previous bin ends just before this one starts.
            if(!bins[i].AsInteger(nextBinStart))
            {
                // bad format, could not turn quality value to integer.
                std::cerr << "Quality Bin Range not an integer: " << bins[i] << std::endl;
                return(-1);
            }
            // Setup this bin.
            binPhredQuals(binStart, nextBinStart - 1);
            // Set the start for the next bin.
            binStart = nextBinStart;
        }
        // Call for the last bin from the last specified start to the MAX_PHRED_QUAL.
        binPhredQuals(nextBinStart, MAX_PHRED_QUAL);
    }
    return(0);
}


int Squeeze::getQualCharFromQemp(uint8_t qemp)
{
    if(!myBinQualS.IsEmpty())
    {
        return(myQualBinMap[qemp + QUAL_CONVERT]);
    }
    // If not squeezing, just return qemp + the conversion.
    return qemp + QUAL_CONVERT;
}


void Squeeze::binPhredQuals(int binStartPhred, int binEndPhred)
{
    // Convert the bin qualities to non-phred.
    int binStart = binStartPhred + QUAL_CONVERT;
    int binEnd = binEndPhred + QUAL_CONVERT;

    if(binEnd > MAX_QUAL_CHAR)
    {
        binEnd = MAX_QUAL_CHAR;
    } 
    // Determine the value for this bin.
    int binValue  = binStart;
    if(myBinMid)
    {
        binValue = (binStart + binEnd)/2;
    }
    else if(myBinHigh)
    {
        binValue = binEnd;
    }

    for(int i = binStart; i <= binEnd; i++)
    {
        myQualBinMap[i] = binValue;
    }
}


void Squeeze::bin(SamRecord& samRecord)
{
    static String qual = "";
    if(!myBinQualS.IsEmpty())
    {
        qual = samRecord.getQuality();
        
        if(qual != "*")
        {
            // Only bin set qualities.
            for(int i = 0; i < qual.Length(); i++)
            {
                // Update the quality by looking its value up
                // in the quality bin map.
                qual[i] = (char)(myQualBinMap[(int)(qual[i])]);
            }
            samRecord.setQuality(qual);
        }
    }
}
