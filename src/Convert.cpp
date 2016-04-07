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

//////////////////////////////////////////////////////////////////////////
// This file contains the processing for the executable option "convert"
// which reads an SAM/BAM file and writes a SAM/BAM file (it can convert 
// between SAM and BAM formats).

#include "Convert.h"
#include "SamFile.h"
#include "Parameters.h"
#include "BgzfFileType.h"
#include "SamValidation.h"

void Convert::printConvertDescription(std::ostream& os)
{
    os << " convert - Convert SAM/BAM to SAM/BAM" << std::endl;
}


void Convert::printDescription(std::ostream& os)
{
    printConvertDescription(os);
}


void Convert::printUsage(std::ostream& os)
{
    BamExecutable::printUsage(os);
    os << "\t./bam convert --in <inputFile> --out <outputFile.sam/bam/ubam (ubam is uncompressed bam)> [--refFile <reference filename>] [--useBases|--useEquals|--useOrigSeq] [--lshift] [--noeof] [--params]" << std::endl;
    os << "\tRequired Parameters:" << std::endl;
    os << "\t\t--in         : the SAM/BAM file to be read" << std::endl;
    os << "\t\t--out        : the SAM/BAM file to be written" << std::endl;
    os << "\tOptional Parameters:" << std::endl;
    os << "\t\t--refFile    : reference file name" << std::endl;
    os << "\t\t--lshift     : left shift indels when writing records\n";
    os << "\t\t--noeof      : do not expect an EOF block on a bam file" << std::endl;
    os << "\t\t--params     : print the parameter settings" << std::endl;
    os << "\t\t--recover    : attempt error recovery while reading a bam file" << std::endl;
    os << "\tOptional Sequence Parameters (only specify one):" << std::endl;
    os << "\t\t--useOrigSeq : Leave the sequence as is (default & used if reference is not specified)" << std::endl;
    os << "\t\t--useBases   : Convert any '=' in the sequence to the appropriate base using the reference (requires --refFile)" << std::endl;
    os << "\t\t--useEquals  : Convert any bases that match the reference to '=' (requires --refFile)" << std::endl;
}

#define SIGNATURE_LENGTH (sizeof(bamRecordStruct))
static bool checkSignature(void *data)
{
    bamRecordStruct *record = (bamRecordStruct *) data;

    if(record->myBlockSize > 60000) return false;

    if(record->myReferenceID < -1 || record->myReferenceID > 1000) return false;

    if(record->myCigarLength > 100) return false;

    if(record->myReadLength > 1000) return false;

    if(record->myMateReferenceID < -1 || record->myMateReferenceID > 1000) return false;

    if(record->myMatePosition < -1) return false;

    if(record->myInsertSize > 100000) return false;

    return true;
}



int Convert::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    String outFile = "";
    String refFile = "";
    bool lshift = false;
    bool noeof = false;
    bool params = false;

    bool useBases = false;
    bool useEquals = false;
    bool useOrigSeq = false;

    bool recover = false;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_STRINGPARAMETER("out", &outFile)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_PARAMETER("lshift", &lshift)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("recover", &recover)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("SequenceConversion")
            EXCLUSIVE_PARAMETER("useBases", &useBases)
            EXCLUSIVE_PARAMETER("useEquals", &useEquals)
            EXCLUSIVE_PARAMETER("useOrigSeq", &useOrigSeq)
        LONG_PHONEHOME(VERSION)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));
    
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
        printUsage(std::cerr);
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--in is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    if(outFile == "")
    {
        printUsage(std::cerr);
        inputParameters.Status();
        // In file was not specified but it is mandatory.
        std::cerr << "--out is a mandatory argument, "
                  << "but was not specified" << std::endl;
        return(-1);
    }

    // Check to see if the ref file was specified.
    // Open the reference.
    GenomeSequence* refPtr = NULL;
    if(refFile != "")
    {
        refPtr = new GenomeSequence(refFile);
    }

    SamRecord::SequenceTranslation translation;
    if((useBases) && (refPtr != NULL))
    {
        translation = SamRecord::BASES;
    }
    else if((useEquals) && (refPtr != NULL))
    {
        translation = SamRecord::EQUAL;
    }
    else
    {
        useOrigSeq = true;
        translation = SamRecord::NONE;
    }
    
    if(params)
    {
        inputParameters.Status();
    }

    // Open the input file for reading.
    SamFile samIn;
    if(recover) samIn.setAttemptRecovery(true);
    samIn.OpenForRead(inFile);

    // Open the output file for writing.
    SamFile samOut;
    samOut.OpenForWrite(outFile);
    samOut.SetWriteSequenceTranslation(translation);
    samOut.SetReference(refPtr);

    // Read the sam header.
    SamFileHeader samHeader;
    samIn.ReadHeader(samHeader);

    // Write the sam header.
    samOut.WriteHeader(samHeader);

    SamRecord samRecord;

    // Set returnStatus to success.  It will be changed
    // to the failure reason if any of the writes fail.
    SamStatus::Status returnStatus = SamStatus::SUCCESS;

    while(1) {
        try {
            // Keep reading records until ReadRecord returns false.
            while(samIn.ReadRecord(samHeader, samRecord))
            {
                // left shift if necessary.
                if(lshift)
                {
                    samRecord.shiftIndelsLeft();
                }

                // Successfully read a record from the file, so write it.
                if(!samOut.WriteRecord(samHeader, samRecord))
                {
                    // Failed to write a record.
                    fprintf(stderr, "%s\n", samOut.GetStatusMessage());
                    returnStatus = samOut.GetStatus();
                }
            }
            break;
        } catch (std::runtime_error e) {
            std::cerr << "Caught runtime error: " << e.what() << "\n";
            if(!recover) {
                std::cerr << "Corrupted BAM file detected - consider using --recover option.\n";
                break;
            }
            std::cerr << "Attempting to resync at next good BGZF block and BAM record.\n";
            // XXX need to resync SamFile stream here
            bool rc = samIn.attemptRecoverySync(checkSignature, SIGNATURE_LENGTH);
            if(rc) {
                std::cerr << "Successful resync - some data lost.\n";
                continue;    // succeeded
            }
            std::cerr << "Failed to re-sync on data stream.\n";
            break;              // failed to resync
        }
    }

    std::cerr << std::endl << "Number of records read = " << 
        samIn.GetCurrentRecordCount() << std::endl;
    std::cerr << "Number of records written = " << 
        samOut.GetCurrentRecordCount() << std::endl;

    if(refPtr != NULL)
    {
        delete(refPtr);
    }

    // Since the reads were successful, return the status based
    // on the status of the writes.  If any failed, return
    // their failure status.
    return(returnStatus);
}


