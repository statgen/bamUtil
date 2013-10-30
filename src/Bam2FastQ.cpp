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

#include "Bam2FastQ.h"
#include "BgzfFileType.h"
#include "BaseUtilities.h"
#include "SamFlag.h"
#include "SamFile.h"
#include "SamHelper.h"


const char* Bam2FastQ::DEFAULT_FIRST_EXT = "/1";
const char* Bam2FastQ::DEFAULT_SECOND_EXT = "/2";



Bam2FastQ::Bam2FastQ()
    : BamExecutable(),
      myPool(),
      myMateMap(true),
      myUnpairedFile(NULL),
      myFirstFile(NULL),
      mySecondFile(NULL),
      myNumMateFailures(0),
      myNumPairs(0),
      myNumUnpaired(0),
      myReverseComp(true),
      myRNPlus(false),
      myFirstRNExt(DEFAULT_FIRST_EXT),
      mySecondRNExt(DEFAULT_SECOND_EXT)
{
}


Bam2FastQ::~Bam2FastQ()
{
    closeFiles();
}


void Bam2FastQ::bam2FastQDescription()
{
    std::cerr << " bam2FastQ - Convert the specified BAM file to fastQs." << std::endl;
}


void Bam2FastQ::description()
{
    bam2FastQDescription();
}


void Bam2FastQ::usage()
{
    BamExecutable::usage();
    std::cerr << "\t./bam bam2FastQ --in <inputFile> [--readName] [--refFile <referenceFile>] [--outBase <outputFileBase>] [--firstOut <1stReadInPairOutFile>] [--secondOut <2ndReadInPairOutFile>] [--unpairedOut <unpairedOutFile>] [--firstRNExt <firstInPairReadNameExt>] [--secondRNExt <secondInPairReadNameExt>] [--rnPlus] [--noReverseComp] [--noeof] [--params]" << std::endl;
    std::cerr << "\tRequired Parameters:" << std::endl;
    std::cerr << "\t\t--in       : the SAM/BAM file to convert to FastQ" << std::endl;
    std::cerr << "\tOptional Parameters:" << std::endl;
    std::cerr << "\t\t--readname      : Process the BAM as readName sorted instead\n"
              << "\t\t                  of coordinate if the header does not indicate a sort order." << std::endl;
    std::cerr << "\t\t--merge         : Generate 1 interleaved (merged) FASTQ for paired-ends (unpaired in a separate file)\n"
              << "\t\t                  use firstOut to override the filename of the interleaved file." << std::endl;
    std::cerr << "\t\t--refFile       : Reference file for converting '=' in the sequence to the actual base" << std::endl;
    std::cerr << "\t\t                  if '=' are found and the refFile is not specified, 'N' is written to the FASTQ" << std::endl;
    std::cerr << "\t\t--firstRNExt    : read name extension to use for first read in a pair\n" 
              << "\t\t                  default is \"" << DEFAULT_FIRST_EXT << "\"\n";
    std::cerr << "\t\t--secondRNExt   : read name extension to use for second read in a pair\n" 
              << "\t\t                  default is \"" << DEFAULT_SECOND_EXT << "\"\n";
    std::cerr << "\t\t--rnPlus        : Add the Read Name/extension to the '+' line of the fastq records\n";
    std::cerr << "\t\t--noReverseComp : Do not reverse complement reads marked as reverse\n";
    std::cerr << "\t\t--noeof         : Do not expect an EOF block on a bam file." << std::endl;
    std::cerr << "\t\t--params        : Print the parameter settings to stderr" << std::endl;
    std::cerr << "\tOptional OutputFile Names:" << std::endl;
    std::cerr << "\t\t--outBase       : Base output name for generated output files" << std::endl;
    std::cerr << "\t\t--firstOut      : Output name for the first in pair file" << std::endl;
    std::cerr << "\t\t                  over-rides setting of outBase" << std::endl;
    std::cerr << "\t\t--secondOut     : Output name for the second in pair file" << std::endl;
    std::cerr << "\t\t                  over-rides setting of outBase" << std::endl;
    std::cerr << "\t\t--unpairedOut   : Output name for unpaired reads" << std::endl;
    std::cerr << "\t\t                  over-rides setting of outBase" << std::endl;
    std::cerr << std::endl;
}


int Bam2FastQ::execute(int argc, char **argv)
{
    // Extract command line arguments.
    String inFile = "";
    bool readName = false;
    String refFile = "";
    String outBase = "";
    String firstOut = "";
    String secondOut = "";
    String unpairedOut = "";

    bool interleave = false;
    bool noeof = false;
    bool params = false;

    myNumMateFailures = 0;
    myNumPairs = 0;
    myNumUnpaired = 0;
    myReverseComp = true;
    myRNPlus = false;
    myFirstRNExt = DEFAULT_FIRST_EXT;
    mySecondRNExt = DEFAULT_SECOND_EXT;

    ParameterList inputParameters;
    BEGIN_LONG_PARAMETERS(longParameterList)
        LONG_PARAMETER_GROUP("Required Parameters")
        LONG_STRINGPARAMETER("in", &inFile)
        LONG_PARAMETER_GROUP("Optional Parameters")
        LONG_PARAMETER("readName", &readName)
        LONG_PARAMETER("merge", &interleave)
        LONG_STRINGPARAMETER("refFile", &refFile)
        LONG_STRINGPARAMETER("firstRNExt", &myFirstRNExt)
        LONG_STRINGPARAMETER("secondRNExt", &mySecondRNExt)
        LONG_PARAMETER("rnPlus", &myRNPlus)
        LONG_PARAMETER("noReverseComp", &myReverseComp)
        LONG_PARAMETER("noeof", &noeof)
        LONG_PARAMETER("params", &params)
        LONG_PARAMETER_GROUP("Optional OutputFile Names")
        LONG_STRINGPARAMETER("outBase", &outBase)
        LONG_STRINGPARAMETER("firstOut", &firstOut)
        LONG_STRINGPARAMETER("secondOut", &secondOut)
        LONG_STRINGPARAMETER("unpairedOut", &unpairedOut)
        END_LONG_PARAMETERS();
   
    inputParameters.Add(new LongParameters ("Input Parameters", 
                                            longParameterList));

    inputParameters.Read(argc-1, &(argv[1]));

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

    // Cannot specify both interleaved & secondOut since secondOut would be N/A.
    if(interleave && !secondOut.IsEmpty())
    {
        usage();
        inputParameters.Status();
        std::cerr << "ERROR: Cannot specify --merge & --secondOut.\n";
        return(-1);
    }


    // Check to see if the out file was specified, if not, generate it from
    // the input filename.
    if(outBase == "")
    {
        // Just remove the extension from the input filename.
        int extStart = inFile.FastFindLastChar('.');
        if(extStart <= 0)
        {
            outBase = inFile;
        }
        else
        {
            outBase = inFile.Left(extStart);
        }
    }

    // Check to see if the first/second/single-ended were specified and
    // if not, set them.
    std::string firstExt = "_1.fastq";
    if(interleave)
    {
        firstExt = "_interleaved.fastq";
    }
    getFileName(firstOut, outBase, firstExt.c_str());
    getFileName(secondOut, outBase, "_2.fastq");
    getFileName(unpairedOut, outBase, ".fastq");

    if(params)
    {
        inputParameters.Status();
    }

    // Open the files for reading/writing.
    // Open prior to opening the output files,
    // so if there is an error, the outputs don't get created.
    SamFile samIn;
    SamFileHeader samHeader;
    samIn.OpenForRead(inFile, &samHeader);

    // Open the output files.
    myUnpairedFile = ifopen(unpairedOut, "w");

    // Only open the first file if it is different than an already opened file.
    if(firstOut != unpairedOut)
    {
        myFirstFile = ifopen(firstOut, "w");
    }
    else
    {
        myFirstFile = myUnpairedFile;
    }

    // If it is interleaved or the 2nd file is not a new name, set it appropriately.
    if(interleave || secondOut == firstOut)
    {
        mySecondFile = myFirstFile;
    }
    else if(secondOut == unpairedOut)
    {
        mySecondFile = myUnpairedFile;
    }
    else
    {
        mySecondFile = ifopen(secondOut, "w");
    }
    
    if(myUnpairedFile == NULL)
    {
        std::cerr << "Failed to open " << unpairedOut
                  << " so can't convert bam2FastQ.\n";
        return(-1);
    }
    if(myFirstFile == NULL)
    {
        std::cerr << "Failed to open " << firstOut
                  << " so can't convert bam2FastQ.\n";
        return(-1);
    }
    if(mySecondFile == NULL)
    {
        std::cerr << "Failed to open " << secondOut
                  << " so can't convert bam2FastQ.\n";
        return(-1);
    }

    if((readName) || (strcmp(samHeader.getSortOrder(), "queryname") == 0))
    {
        readName = true;
    }
    else
    {
        // defaulting to coordinate sorted.
        samIn.setSortedValidation(SamFile::COORDINATE);
    }

    // Setup the '=' translation if the reference was specified.
    if(!refFile.IsEmpty())
    {
        GenomeSequence* refPtr = new GenomeSequence(refFile);
        samIn.SetReadSequenceTranslation(SamRecord::BASES);
        samIn.SetReference(refPtr);
    }

    SamRecord* recordPtr;
    int16_t samFlag;

    SamStatus::Status returnStatus = SamStatus::SUCCESS;
    while(returnStatus == SamStatus::SUCCESS)
    {
        recordPtr = myPool.getRecord();
        if(recordPtr == NULL)
        {
            // Failed to allocate a new record.
            throw(std::runtime_error("Failed to allocate a new SAM/BAM record"));
        }
        if(!samIn.ReadRecord(samHeader, *recordPtr))
        {
            // Failed to read a record.
            returnStatus = samIn.GetStatus();
            continue;
        }

        // Have a record.  Check to see if it is a pair or unpaired read.
        samFlag = recordPtr->getFlag();
        if(SamFlag::isPaired(samFlag))
        {
            if(readName)
            {
                handlePairedRN(*recordPtr);
            }
            else
            {
                handlePairedCoord(*recordPtr);
            }
        }
        else
        {
            ++myNumUnpaired;
            writeFastQ(*recordPtr, myUnpairedFile);
        }
    }

    // Flush All
    cleanUpMateMap(0, true);

    if(returnStatus == SamStatus::NO_MORE_RECS)
    {
        returnStatus = SamStatus::SUCCESS;
    }

    samIn.Close();
    closeFiles();
    
    // Output the results
    std::cerr << "\nFound " << myNumPairs << " read pairs.\n";
    std::cerr << "Found " << myNumUnpaired << " unpaired reads.\n";
    if(myNumMateFailures != 0)
    {
        std::cerr << "Failed to find mates for " << myNumMateFailures
                  << " reads, so they were written as unpaired\n"
                  << "  (not included in either of the above counts).\n";
    }

    return(returnStatus);
}


void Bam2FastQ::handlePairedRN(SamRecord& samRec)
{
    static SamRecord* prevRec = NULL;
    static std::string prevRN = "";

    if(prevRec == NULL)
    {
        prevRec = &samRec;
    }
    else
    {
        if(strcmp(prevRec->getReadName(), samRec.getReadName()) != 0)
        {
            // Read Name does not match, error, did not find pair.
            std::cerr << "Paired Read, " << prevRec->getReadName()
                      << " but couldn't find mate, so writing as "
                      << "unpaired (single-ended)\n";
            ++myNumMateFailures;
            writeFastQ(*prevRec, myUnpairedFile);
            // Save this record to check against the next one.
            prevRec = &samRec;
        }
        else
        {
            // Matching ReadNames.
            // Found the mate.
            ++myNumPairs;
            // Check which is the first in the pair.
            if(SamFlag::isFirstFragment(samRec.getFlag()))
            {
                if(SamFlag::isFirstFragment(prevRec->getFlag()))
                {
                    std::cerr << "Both reads of " << samRec.getReadName()
                              << " are first fragment, so "
                              << "splitting one to be in the 2nd fastq.\n";
                }
                writeFastQ(samRec, myFirstFile, myFirstRNExt.c_str());
                writeFastQ(*prevRec, mySecondFile, mySecondRNExt.c_str());
            }
            else
            {
                if(!SamFlag::isFirstFragment(prevRec->getFlag()))
                {
                    std::cerr << "Neither read of " << samRec.getReadName()
                              << " are first fragment, so "
                              << "splitting one to be in the 2nd fastq.\n";
                }
                writeFastQ(*prevRec, myFirstFile, myFirstRNExt.c_str());
                writeFastQ(samRec, mySecondFile, mySecondRNExt.c_str());
            }
            // No previous record.
            prevRec = NULL;
        }
    }
}


void Bam2FastQ::handlePairedCoord(SamRecord& samRec)
{
    static uint64_t readPos;
    static uint64_t matePos;
    static SamRecord* mateRec;

    // This is a paired record, so check for its mate.
    readPos = SamHelper::combineChromPos(samRec.getReferenceID(),
                                         samRec.get0BasedPosition());
    matePos = SamHelper::combineChromPos(samRec.getMateReferenceID(), 
                                         samRec.get0BasedMatePosition());
 
    // Check to see if the mate is prior to this record.
    if(matePos <= readPos)
    {
        // The mate is prior to this position, so look for this record.
        mateRec = myMateMap.getMate(samRec);
        if(mateRec == NULL)
        {
            // If they are the same position, add it to the map.
            if(matePos == readPos)
            {
                myMateMap.add(samRec);

                // Check to see if the mate map can be cleaned up prior
                // to this position.
                cleanUpMateMap(readPos);
            }
            else
            {
                // Paired Read, but could not find mate.
                std::cerr << "Paired Read, " << samRec.getReadName()
                          << " but couldn't find mate, so writing as "
                          << "unpaired (single-ended)\n";
                ++myNumMateFailures;
                writeFastQ(samRec, myUnpairedFile);
            }
        }
        else
        {
            // Found the mate.
            ++myNumPairs;
            // Check which is the first in the pair.
            if(SamFlag::isFirstFragment(samRec.getFlag()))
            {
                if(SamFlag::isFirstFragment(mateRec->getFlag()))
                {
                    std::cerr << "Both reads of " << samRec.getReadName()
                              << " are first fragment, so "
                              << "splitting one to be in the 2nd fastq.\n";
                }
                writeFastQ(samRec, myFirstFile, myFirstRNExt.c_str());
                writeFastQ(*mateRec, mySecondFile, mySecondRNExt.c_str());
            }
            else
            {
                if(!SamFlag::isFirstFragment(mateRec->getFlag()))
                {
                    std::cerr << "Neither read of " << samRec.getReadName()
                              << " are first fragment, so "
                              << "splitting one to be in the 2nd fastq.\n";
                }
                writeFastQ(*mateRec, myFirstFile, myFirstRNExt.c_str());
                writeFastQ(samRec, mySecondFile, mySecondRNExt.c_str());
            }
        }
    }
    else
    {
        // Haven't gotten to the mate yet, so store it.
        myMateMap.add(samRec);

        // Check to see if the mate map can be cleaned up.
        cleanUpMateMap(readPos);
    }
}


void Bam2FastQ::writeFastQ(SamRecord& samRec, IFILE filePtr,
                             const char* readNameExt)
{
    static int16_t flag;
    static std::string sequence;
    static String quality;

    if(filePtr == NULL)
    {
        return;
    }

    flag = samRec.getFlag();
    const char* readName = samRec.getReadName();
    sequence = samRec.getSequence();
    quality = samRec.getQuality();
    
    if(SamFlag::isReverse(flag) && myReverseComp)
    {
        // It is reverse, so reverse compliment the sequence
        BaseUtilities::reverseComplement(sequence);
        // Reverse the quality.
        quality.Reverse();
    }
    else
    {
        // Ensure it is all capitalized.
        int seqLen = sequence.size();
        for (int i = 0; i < seqLen; i++)
        {
            sequence[i] = (char)toupper(sequence[i]);
        }
    }
    
    if(myRNPlus)
    {

        ifprintf(filePtr, "@%s%s\n%s\n+%s%s\n%s\n", readName, readNameExt,
                 sequence.c_str(), readName, readNameExt, quality.c_str());
    }
    else
    {
        ifprintf(filePtr, "@%s%s\n%s\n+\n%s\n", readName, readNameExt,
                 sequence.c_str(), quality.c_str());
    }
    // Release the record.
    myPool.releaseRecord(&samRec);
}


void Bam2FastQ::cleanUpMateMap(uint64_t readPos, bool flushAll)
{
    static SamRecord* firstRec;

    firstRec = myMateMap.first();
    while(firstRec != NULL)
    {
        uint64_t firstChromPos = 
            SamHelper::combineChromPos(firstRec->getMateReferenceID(),
                                       firstRec->get0BasedMatePosition());
        if((firstChromPos < readPos) || flushAll)
        {
            // Already past the mate's position, so did not find the mate for
            // this record, write it out as single-ended.
            std::cerr << "Paired Read, " << firstRec->getReadName()
                      << " but couldn't find mate, so writing as "
                      << "unpaired (single-ended)\n";
            writeFastQ(*firstRec, myUnpairedFile);
            ++myNumMateFailures;

            // Remove this record from the map.
            myMateMap.popFirst();
            firstRec = myMateMap.first();
        }
        else
        {
            // We have not yet hit the mate positions in the map, so
            // no more records to release.
            firstRec = NULL;
        }
    }
}


void Bam2FastQ::closeFiles()
{
    // NULL out any duplicate file pointers
    // so files are only closed once.
    if(myFirstFile == myUnpairedFile)
    {
        myFirstFile = NULL;
    }
    if(mySecondFile == myUnpairedFile)
    {
        mySecondFile = NULL;
    }
    if(mySecondFile == myFirstFile)
    {
        mySecondFile = NULL;
    }

    if(myUnpairedFile != NULL)
    {
        ifclose(myUnpairedFile);
        myUnpairedFile = NULL;
    }
    if(myFirstFile != NULL)
    {
        ifclose(myFirstFile);
        myFirstFile = NULL;
    }
    if(mySecondFile != NULL)
    {
        ifclose(mySecondFile);
        mySecondFile = NULL;
    }
}


void Bam2FastQ::getFileName(String& fn, const String& outBase, const char* ext)
{
    if(fn.IsEmpty())
    {
        if(outBase[0] != '-')
        {
            fn = outBase + ext;
        }
        else
        {
            // starts with a '-', so writing to stdout, don't change the name.
            fn = outBase;
        }
    }
}


