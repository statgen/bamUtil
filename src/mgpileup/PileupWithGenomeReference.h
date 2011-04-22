/*
 *  Copyright (C) 2010-2011  Regents of the University of Michigan
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

#ifndef __PILEUP_WITH_GENOME_REFERENCE_H__
#define __PILEUP_WITH_GENOME_REFERENCE_H__

#include "Pileup.h"
#include "GenomeSequence.h"
#include "InputFile.h"
#include "PreviousPileup.h"

struct Region
{
    std::string chrom;
    uint32_t start;
    uint32_t end;
    vector<int> positions;
    unsigned int currentPosition;
};

template <class PILEUP_TYPE, 
          class FUNC_CLASS = defaultPileup<PILEUP_TYPE> >
class PileupWithGenomeReference:public Pileup<PILEUP_TYPE, FUNC_CLASS>{

public:
    PileupWithGenomeReference(const std::string& refSeqFileName,
                              bool addDelAsBase = false,
                              const FUNC_CLASS& fp = FUNC_CLASS());

    PileupWithGenomeReference(int window, const std::string& refSeqFileName,
                              bool addDelAsBase = false,
                              const FUNC_CLASS& fp = FUNC_CLASS());
           
    virtual int processFile(const std::string& bamFileName,  
                            const std::string& outputVCFFileName,
                            uint16_t excludeFlag = 0x0704, 
                            uint16_t includeFlag = 0);        

    virtual int processFile(const std::string& bamFileName, 
                            const std::string& inputVCFFileName,
                            const std::string& prevPileupFileName,
                            const std::string& outputVCFFileName,
                            uint32_t maxStoredPrevPileupLines,
                            uint16_t excludeFlag = 0x0704,
                            uint16_t includeFlag = 0);
                    
    virtual void processAlignment(SamRecord& record, Region& region);

    virtual void analyzeHead();

    void initLogGLMatrix();
    
private:
    int processRegion(Region& currentRegion, 
                      uint16_t excludeFlag,
                      uint16_t includeFlag);
    void resetElement(PILEUP_TYPE& element, int position);


    bool myAddDelAsBase;
    GenomeSequence myRefSeq;
    VcfFile myOutputVCFFile;

    SamFile mySamIn;
    SamFileHeader mySamHeader;

    PreviousPileup myPrevPileup;

    double*** myLogGLMatrix;
};

template <class PILEUP_TYPE, class FUNC_CLASS>
PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::PileupWithGenomeReference(const std::string& refSeqFileName,
                                                                              bool addDelAsBase, 
                                                                              const FUNC_CLASS& fp)
    : 	Pileup<PILEUP_TYPE, FUNC_CLASS>(),
        myAddDelAsBase(addDelAsBase),
        myRefSeq(refSeqFileName.c_str(), false),
        myOutputVCFFile(),
        mySamIn(),
        mySamHeader(),
        myPrevPileup()
{
    initLogGLMatrix();
}

template <class PILEUP_TYPE, class FUNC_CLASS>
PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::PileupWithGenomeReference(int window, 
                                                                              const std::string& refSeqFileName,
                                                                              bool addDelAsBase,
                                                                              const FUNC_CLASS& fp)
    :	Pileup<PILEUP_TYPE, FUNC_CLASS>(window),
        myAddDelAsBase(addDelAsBase),
        myRefSeq(refSeqFileName.c_str()),
        myOutputVCFFile(),
        mySamIn(),
        mySamHeader(),
        myPrevPileup()
{
    initLogGLMatrix();
}

template <class PILEUP_TYPE, class FUNC_CLASS>
int PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processFile(const std::string& bamFileName,
                                                                    const std::string& outputVCFFileName, 
                                                                    uint16_t excludeFlag,
                                                                    uint16_t includeFlag)
{
    /////////////////////////////////////
    // Open the output files.
    if(!myOutputVCFFile.openForWrite(outputVCFFileName.c_str()))
    {
        abort();
    }

    Pileup<PILEUP_TYPE, FUNC_CLASS>::processFile(bamFileName, excludeFlag, includeFlag);

    myOutputVCFFile.close();

    return(0);  
}

template <class PILEUP_TYPE, class FUNC_CLASS>
int PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processFile(const std::string& bamFileName, 
                                                                    const std::string& inputVCFFileName, 
                                                                    const std::string& prevPileupFileName,
                                                                    const std::string& outputVCFFileName, 
                                                                    uint32_t maxStoredPrevPileupLines,
                                                                    uint16_t excludeFlag,
                                                                    uint16_t includeFlag)
{
    /////////////////////////////////////
    // Open the input and output files.
    if(!myOutputVCFFile.openForWrite(outputVCFFileName.c_str()))
    {
        abort();
    }

    VcfFile vcfRegions;
    if(!vcfRegions.openForRead(inputVCFFileName.c_str()))
    {
        abort();
    }

    myPrevPileup.open(prevPileupFileName, maxStoredPrevPileupLines);

    // Open the bam file.
    bool hasBamFile = false;

    if(bamFileName.length() != 0)
    {
        // There is a bam file, so open it.
        if(!mySamIn.OpenForRead(bamFileName.c_str()))
        {
            fprintf(stderr, "%s\n", mySamIn.GetStatusMessage());
        }
        else if(!mySamIn.ReadHeader(mySamHeader) || !mySamIn.ReadBamIndex())
        {
            fprintf(stderr, "%s\n", mySamIn.GetStatusMessage());
            // Couldn't read, so close the file.
            mySamIn.Close();
        }
        else
        {
            // The file needs to be sorted by coordinate.
            mySamIn.setSortedValidation(SamFile::COORDINATE);
            hasBamFile = true;
        }
    }

    //////////////////////////////////////
    // Now that the files are open, process them.
    Region currentRegion;
    currentRegion.chrom = "";
    currentRegion.start = 0;
    currentRegion.end = 0;
    currentRegion.positions.clear();
    currentRegion.currentPosition = 0;

    std::string currentChrom = "";
    int currentChromID = -1;

    // Read until there are no more chromosome/position lines in the input vcf position file.
    VcfFile::VcfDataLine nextPos;

    //////////////////////////////////////////////////////////////
    // Loop through getting the positions that should be piled-up.
    int status = 0;
    while(vcfRegions.getNextDataLine(nextPos) && (status == 0))
    {
        // If a region is set, check to see if it is in the same region.
        // It is not in the same region if:
        //   a) it is on a different chromosome.
        //   b) it is 1000 bases or more past the end of this region.
        if((currentRegion.positions.size() != 0) &&
           ((nextPos.chromosome!=currentRegion.chrom) ||
            (nextPos.position-currentRegion.end>=1000)))
        {
            // Not in the same region, so process the region.
            status = processRegion(currentRegion, excludeFlag, includeFlag);
        }
        
        if(nextPos.chromosome != currentChrom)
        {
            currentChrom = nextPos.chromosome;
            currentChromID = mySamHeader.getReferenceID(currentChrom.c_str());
        }
        
        // Determine where the data for this position is to come from.
        // Check it against the Previous VCF.
        if(myPrevPileup.findTarget(currentChromID, nextPos.position, mySamHeader))
        {
            // Found this position in the previous piluep.
            if(currentRegion.positions.size() == 0)
            {
                // No bam entries yet, so just write this entry.
                myPrevPileup.writeCurrentLine(myOutputVCFFile);
            }
            else
            {
                // Store the line for later.
                if(!myPrevPileup.storeCurrent())
                {
                    // The previous pileup has reached its max storage size, 
                    // so process the region to write those lines.
                    status = processRegion(currentRegion, excludeFlag, includeFlag);
                    if(status == 0)
                    {
                        // Just wrote any bam entries, so can just write this line.
                        myPrevPileup.writeCurrentLine(myOutputVCFFile);
                    }
                }
            }
        }
        else
        {
            // Need to read from the bam file.
            if(currentRegion.positions.size() == 0)
            {
                //create new current region
                currentRegion.chrom=nextPos.chromosome;
                currentRegion.start=nextPos.position;
                currentRegion.currentPosition = 0;
            }
            //extend region
            currentRegion.end = nextPos.position;
            currentRegion.positions.push_back(nextPos.position);
        }
    }

    // If a region was set, process the last region.
    if(currentRegion.positions.size() != 0)
    {
        // process the last region
        status = processRegion(currentRegion, excludeFlag, includeFlag);
    }
    
    myOutputVCFFile.close();
    vcfRegions.close();

    return(status);
}


template <class PILEUP_TYPE, class FUNC_CLASS>
int PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processRegion(Region& currentRegion, 
                                                                      uint16_t excludeFlag,
                                                                      uint16_t includeFlag)
{
    // Set the read section.
    if(!mySamIn.SetReadSection(mySamHeader.getReferenceID(currentRegion.chrom.c_str()), 
                             currentRegion.start-1, currentRegion.end))
    {
        std::cerr << "chrom:" << currentRegion.chrom << ":" <<  currentRegion.start 
                  << "-" << currentRegion.end << "\n";
    }
    else
    {
        // Iterate over all records
        SamRecord record;
        while (mySamIn.ReadRecord(mySamHeader, record))
        {
            uint16_t flag = record.getFlag();
            if(flag & excludeFlag)
            {
                // This record has an excluded flag set, 
                // so continue to the next one.
                continue;
            }
            if((flag & includeFlag) != includeFlag)
            {
                // This record does not have all required flags set, 
                // so continue to the next one.
                continue;
            }
            
            processAlignment(record, currentRegion);
        }
    }

    Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup();
    
    // Flush the rest of the vcf queue.
    myPrevPileup.writeStored(myOutputVCFFile);

    int returnValue = 0;
    if(mySamIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", mySamIn.GetStatusMessage());
        returnValue = mySamIn.GetStatus();
    }

    currentRegion.positions.clear();

    return(returnValue);  
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processAlignment(SamRecord& record, Region& region)
{
    int refPosition = record.get0BasedPosition();
    int refID = record.getReferenceID();

    // Flush any elements from the pileup that are prior to this record
    // since the file is sorted, we are done with those positions.
    Pileup<PILEUP_TYPE, FUNC_CLASS>::flushPileup(refID, refPosition);
    
    //search for first location in region.positions that is >= the start position of the record
    while (region.positions.at(region.currentPosition)-1 < refPosition)
    {
        // Increment the current position because all future records are beyond this position.
        ++(region.currentPosition);
    }

    // For each position til the end of the positions or end of the alignment, add this record.
    for (unsigned int k = region.currentPosition; 
         ((k < region.positions.size()) &&
          (region.positions.at(k)-1 <= record.get0BasedAlignmentEnd()));
         ++k)
    {
        Pileup<PILEUP_TYPE, FUNC_CLASS>::addAlignmentPosition(region.positions.at(k)-1, record);
    }
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::analyzeHead()
{
    myPrevPileup.writeStoredUpToPos(Pileup<PILEUP_TYPE, FUNC_CLASS>::pileupHead + 1, myOutputVCFFile);
    
    Pileup<PILEUP_TYPE, FUNC_CLASS>::analyzeHead();
}


template <class PILEUP_TYPE, class FUNC_CLASS>
void PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::resetElement(PILEUP_TYPE& element,
                                                                      int position)
{
    element.reset(position, &myRefSeq, myOutputVCFFile, myAddDelAsBase, myLogGLMatrix);
}

template <class PILEUP_TYPE, class FUNC_CLASS>
void PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::initLogGLMatrix()
{
    myLogGLMatrix = (double***)malloc(50*sizeof(double**));
    for (int i=0; i<50; ++i)
    {
        myLogGLMatrix[i] = (double**)malloc(10*sizeof(double*));
        for (int j=0; j<10; ++j)
        {
            myLogGLMatrix[i][j] = (double*)malloc(4*sizeof(double));
        }
    }
 
    std::string genotypes[10] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
    char bases[4] = {'A','C','G','T'};
    double e[50];
    char allele1;
    char allele2;
    char base;		

    for (int i=0; i<=49; ++i)
    {	
        e[i] = pow(10, -i/10.0)/3;

        for (int j=0; j<=9; ++j)
    	{
            for (int k=0; k<=3; ++k)
            {  	
                allele1 = genotypes[j].c_str()[0];
                allele2 = genotypes[j].c_str()[1];
                base = bases[k];

                if(allele1==allele2)
                {
                    myLogGLMatrix[i][j][k] = (base==allele1) ? log10(1-e[i]*3) : log10(e[i]);
                }
                else
                {
                    myLogGLMatrix[i][j][k] = (base==allele1 || base==allele2) ? log10(0.5-e[i]) : log10(e[i]);
                }
            }
    	}
    }
}
#endif
