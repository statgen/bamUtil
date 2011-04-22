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

#include <stdexcept>
#include <stdint.h>
#include <stdarg.h>
#include <math.h>

#include "PileupElementBaseQual.h"
#include "GenomeSequence.h"
#include "BaseUtilities.h"

PileupElementBaseQual::PileupElementBaseQual()
    : PileupElement(),
      myBases(NULL),
      myMapQualities(NULL),
      myQualities(NULL),
      myStrands(NULL),
      myCycles(NULL),
      myGLScores(NULL),
      myAllocatedSize(0),
      myIndex(-1),
      myAddDelAsBase(false),
      myRefSeq(NULL),
      myVcfOutFile(NULL)
{
    myAllocatedSize = 1024;
    myBases = (char*)malloc(myAllocatedSize + 1);
    myMapQualities = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myQualities = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myStrands = (char*)malloc(myAllocatedSize + 1);
    myCycles = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myGLScores = (int16_t*)malloc(myAllocatedSize * sizeof(int16_t));
    if((myBases == NULL ) 
       || (myMapQualities == NULL)
       || (myQualities == NULL)
       || (myStrands == NULL)
       || (myCycles == NULL)
       || (myGLScores == NULL))
    {                     
        // TODO, check for malloc failures.
        std::cerr << "Failed Memory Allocation\n";
    }
}                         

PileupElementBaseQual::PileupElementBaseQual(bool addDelAsBase)
    : PileupElement(),
      myBases(NULL),
      myMapQualities(NULL),
      myQualities(NULL),
      myStrands(NULL),
      myCycles(NULL),
      myGLScores(NULL),
      myAllocatedSize(0),
      myIndex(-1),
      myAddDelAsBase(addDelAsBase),
      myRefSeq(NULL),
      myVcfOutFile(NULL)
{
    myAllocatedSize = 1024;
    myBases = (char*)malloc(myAllocatedSize + 1);
    myMapQualities = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myQualities = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myGLScores = (int16_t*)malloc(myAllocatedSize * sizeof(int16_t));
    myStrands = (char*)malloc(myAllocatedSize + 1);
    myCycles = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    if((myBases == NULL ) 
       || (myMapQualities == NULL)
       || (myQualities == NULL)
       || (myStrands == NULL)
       || (myCycles == NULL)
       || (myGLScores == NULL))
    {                     
        // TODO, check for malloc failures.
        std::cerr << "Failed Memory Allocation\n";
    }                     
}  
                          
PileupElementBaseQual::PileupElementBaseQual(const PileupElementBaseQual& q)
    : PileupElement(),    
      myBases(NULL),      
      myMapQualities(NULL),
      myQualities(NULL), 
      myStrands(NULL),    
      myCycles(NULL),    
      myGLScores(NULL), 
      myAllocatedSize(0), 
      myIndex(-1),
      myRefSeq(NULL),
      myVcfOutFile(NULL)
{                         
    myAllocatedSize = 1024;
    myBases = (char*)malloc(myAllocatedSize + 1);
    myMapQualities = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myQualities = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myStrands = (char*)malloc(myAllocatedSize + 1);
    myCycles = (int8_t*)malloc(myAllocatedSize * sizeof(int8_t));
    myGLScores = (int16_t*)malloc(myAllocatedSize * sizeof(int16_t));
    myAddDelAsBase = q.myAddDelAsBase;
    if((myBases == NULL ) 
       || (myMapQualities == NULL)
       || (myQualities == NULL)
       || (myStrands == NULL)
       || (myCycles == NULL)
       || (myGLScores == NULL))
    {
        // TODO, check for malloc failures.
        std::cerr << "Failed Memory Allocation\n";
    }
}

PileupElementBaseQual::~PileupElementBaseQual()
{
    if(myBases != NULL)
    {
        free(myBases);
        myBases = NULL;
    }
    if(myMapQualities != NULL)
    {
        free(myMapQualities);
        myQualities = NULL;
    }
    if(myQualities != NULL)
    {
        free(myQualities);
        myQualities = NULL;
    }
    if(myStrands != NULL)
    {
        free(myStrands);
        myStrands = NULL;
    }
    if(myCycles != NULL)
    {
        free(myCycles);
        myCycles = NULL;
    }
    if(myGLScores != NULL)
    {
        free(myGLScores);
        myGLScores = NULL;
    }
}

void PileupElementBaseQual::computeGLScores(int index, int16_t* GLScores, char* bases, int8_t* baseQualities)
{
    int base;
    double result;
    for (int genotype=0; genotype<=9; ++genotype)
    {
    	result = 0;
        for (int i=0; i<=myIndex; ++i)
    	{
            switch(bases[i])
            {
                case 'A':
                    base = 0;
                    break;
                case 'C':
                    base = 1;
                    break;
                case 'G':
                    base = 2;
                    break;
                case 'T':
                    base = 3;
                    break;    			
                default :
                    base = -1;
            }

            if(base!=-1)
            {
                result += myLogGLMatrix[baseQualities[i]][genotype][base];
            }
    	}

    	GLScores[genotype] = -10*result;
    }
}

const char* PileupElementBaseQual::getRefAllele() 
{ 
    return(myRefAllele.c_str()); 
}
    
// Add an entry to this pileup element.  
void PileupElementBaseQual::addEntry(SamRecord& record)
{
    // Call the base class:
    PileupElement::addEntry(record);

    if(myRefAllele.empty())
    {
    	genomeIndex_t markerIndex = (*myRefSeq).getGenomePosition(getChromosome(), static_cast<uint32_t>(getRefPosition()+1));
        myRefAllele = (*myRefSeq)[markerIndex];
    }

    // Increment the index
    ++myIndex;
    
    // if the index has gone beyond the allocated space, double the size.
    if(myIndex >= myAllocatedSize)
    {
        char* tempBuffer = (char*)realloc(myBases, myAllocatedSize * 2);
        if(tempBuffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myBases = tempBuffer;
        int8_t* tempInt8Buffer = (int8_t*)realloc(myMapQualities, myAllocatedSize * 2 * sizeof(int8_t));
        if(tempInt8Buffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myMapQualities = tempInt8Buffer; 
        tempInt8Buffer = (int8_t*)realloc(myQualities, myAllocatedSize * 2 * sizeof(int8_t));
        if(tempInt8Buffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myQualities = tempInt8Buffer;
        tempBuffer = (char*)realloc(myStrands, myAllocatedSize * 2);
        if(tempBuffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myStrands = tempBuffer;
        tempInt8Buffer = (int8_t*)realloc(myCycles, myAllocatedSize * 2 * sizeof(int8_t));
        if(tempInt8Buffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myCycles = tempInt8Buffer; 
        int16_t* tempInt16Buffer = (int16_t*)realloc(myGLScores, myAllocatedSize * 2 * sizeof(int16_t));
        if(tempInt8Buffer == NULL)
        {
            std::cerr << "Memory Allocation Failure\n";
            // TODO
            return;
        }
        myGLScores = tempInt16Buffer;
        myAllocatedSize = myAllocatedSize * 2;
    }

    Cigar* cigar = record.getCigarInfo();
    
    if(cigar == NULL)
    {
        throw std::runtime_error("Failed to retrieve cigar info from the record.");
    }

    int32_t readIndex = 
        cigar->getQueryIndex(getRefPosition(), record.get0BasedPosition());

    // If the readPosition is N/A, this is a deletion.
    if(readIndex != CigarRoller::INDEX_NA)
    {
        char base = record.getSequence(readIndex);
        int8_t mapQual = record.getMapQuality();
        //-33 to obtain the PHRED base quality
        char qual = record.getQuality(readIndex);
        if(qual == BaseUtilities::UNKNOWN_QUALITY_CHAR)
        {
            qual = ' ';
        }
        else
        {
            qual -= 33;
        }
        char strand = (record.getFlag() & 0x0010) ? 'R' : 'F';
        int cycle = strand == 'F' ? readIndex + 1 : record.getReadLength() -  readIndex;
        myBases[myIndex] = base;
        myMapQualities[myIndex] = mapQual;
        myQualities[myIndex] = qual;
        myStrands[myIndex] = strand;
        myCycles[myIndex] = cycle;
    }
    else if(myAddDelAsBase)
    {
        int8_t mapQual = record.getMapQuality();
        char strand = (record.getFlag() & 0x0010) ? 'R' : 'F';
        myBases[myIndex] = '-';
        myMapQualities[myIndex] = mapQual;
        myQualities[myIndex] = -1;
        myStrands[myIndex] = strand;
        myCycles[myIndex] = -1;
    }
    else
    {
        // Do not add a deletion.
        // Did not add any entries, so decrement the index counter since the
        // index was not used.
        --myIndex;
    }
}

void PileupElementBaseQual::analyze()
{
    if(getRefPosition() != UNSET_POSITION && myIndex != -1)
    {
    	char tempCStr[11];
    	std::string tempStr;
    	tempStr.append(getChromosome());
    	tempStr.append("\t");
    	sprintf(tempCStr, "%d", getRefPosition()+1 );
    	tempStr.append(tempCStr);
    	tempStr.append("\t.\t");
    	tempStr.append(getRefAllele());
    	tempStr.append("\t.\t.\t.\t.\t");
         
        tempStr.append("N:BASE:MAPQ:BASEQ:STRAND:CYCLE:GL\t");
        
        sprintf(tempCStr, "%d", myIndex+1 );
    	tempStr.append(tempCStr);
    	tempStr.append(":");
        
        sprintf(tempCStr, "%c", myBases[0]);
        tempStr.append(tempCStr);
        for (int i=1; i<=myIndex; ++i)
        {
            sprintf(tempCStr, ",%c", myBases[i]);
            tempStr.append(tempCStr);
        }
        tempStr.append(":");	
		
        sprintf(tempCStr, "%d", myMapQualities[0]);
        tempStr.append(tempCStr);
        for (int i=1; i<=myIndex; ++i)
        {
            sprintf(tempCStr, ",%d", myMapQualities[i]);
            tempStr.append(tempCStr);
        }
        tempStr.append(":");

        sprintf(tempCStr, "%d", myQualities[0]);
        tempStr.append(tempCStr);    
        for (int i=1; i<=myIndex; ++i)
        {
            sprintf(tempCStr, ",%d", myQualities[i]);
            tempStr.append(tempCStr);
        }
        tempStr.append(":");
 
        sprintf(tempCStr, "%c", myStrands[0]);
        tempStr.append(tempCStr);
        for (int i=1; i<=myIndex; ++i)
        {
            sprintf(tempCStr, ",%c", myStrands[i]);
            tempStr.append(tempCStr);
        }
        tempStr.append(":");
 
        sprintf(tempCStr, "%d", myCycles[0]);
        tempStr.append(tempCStr); 
        for (int i=1; i<=myIndex; ++i)
        {
            sprintf(tempCStr, ",%d", myCycles[i]);
            tempStr.append(tempCStr);        	
        }
        tempStr.append(":");
         
        computeGLScores(myIndex, myGLScores, myBases, myQualities);
        
        sprintf(tempCStr, "%d", myGLScores[0]);
        tempStr.append(tempCStr);    
        for (int i=1; i<=9; ++i)
        {
            sprintf(tempCStr, ",%d", myGLScores[i]);
            tempStr.append(tempCStr);
        }
        
        tempStr.append("\n");
       	myVcfOutFile->writeLine(tempStr);
    }
    
    //to ensure this does not print when reflushed
    myIndex=-1;
}


void PileupElementBaseQual::reset(int refPosition)
{	
    // Call the base class.      
    PileupElement::reset(refPosition);

    myIndex = -1;
}


void PileupElementBaseQual::reset(int refPosition, GenomeSequence* refSeq, VcfFile& vcfOutFile, bool addDelAsBase, double*** logGLMatrix)
{
    // Assign pointer to myLogGLMatrix
    if (myLogGLMatrix == NULL)
    {
        myLogGLMatrix = logGLMatrix;
    }
		
    // Assign pointer to myRefSeq
    if (myRefSeq == NULL)
    {
        myRefSeq = refSeq;
    }

    // Assign pointer to myVcfOutFile
    if (myVcfOutFile == NULL)
    {
        myVcfOutFile = &vcfOutFile;
    }
	
    myAddDelAsBase = addDelAsBase;	
    myRefAllele.clear();
	
    // Call the base class.      
    PileupElement::reset(refPosition);

    myIndex = -1;
}

