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

#ifndef __PILEUP_WITH_GENOME_REFERENCE_H__
#define __PILEUP_WITH_GENOME_REFERENCE_H__

#include "Pileup.h"
#include "GenomeSequence.h"
#include "InputFile.h"

struct Region
{
	std::string chrom;
	uint32_t start;
	uint32_t end;
	vector<uint32_t> *positions;
	uint32_t currentPosition;
};

template <class PILEUP_TYPE, 
          class FUNC_CLASS = defaultPileup<PILEUP_TYPE> >
class PileupWithGenomeReference:public Pileup<PILEUP_TYPE, FUNC_CLASS>{

public:
	PileupWithGenomeReference(const std::string& refSeqFileName,
							  bool addDelAsBase = false, 
							  bool inputVCFFileIsGZipped = false,
							  bool outputVCFFileIsGZipped = false, 
		                      const FUNC_CLASS& fp = FUNC_CLASS());

    PileupWithGenomeReference(int window, const std::string& refSeqFileName,
    	   					  bool addDelAsBase = false,  
							  bool inputVCFFileIsGZipped = false,
							  bool outputVCFFileIsGZipped = false, 
    					      const FUNC_CLASS& fp = FUNC_CLASS());
           
    virtual int processFile(const std::string& bamFileName,  
                    const std::string& outputVCFFileName,
                    uint16_t excludeFlag = 0x0704, 
                    uint16_t includeFlag = 0);        

	virtual int processFile(const std::string& bamFileName, 
					const std::string& bamIndexFileName, 
					const std::string& inputVCFFileName,
                    const std::string& outputVCFFileName,
                    uint16_t excludeFlag = 0x0704,
                    uint16_t includeFlag = 0);
                    
    virtual void processAlignment(SamRecord& record, Region* region);
    	                                             					                    
    void initLogGLMatrix();
    
private:
	void resetElement(PILEUP_TYPE& element, int position);
	
	bool myAddDelAsBase;
	bool inputVCFFileIsGZipped;
	bool outputVCFFileIsGZipped;
	GenomeSequence myRefSeq;
	InputFile* myOutputVCFFile;
	double*** myLogGLMatrix;
};

template <class PILEUP_TYPE, class FUNC_CLASS>
PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::PileupWithGenomeReference(const std::string& refSeqFileName,
																		      bool addDelAsBase, 
							  											      bool inputVCFFileIsGZipped,
							  												  bool outputVCFFileIsGZipped, 
																		      const FUNC_CLASS& fp)
: 	Pileup<PILEUP_TYPE>(),
	myAddDelAsBase(addDelAsBase),
	inputVCFFileIsGZipped(inputVCFFileIsGZipped),
	outputVCFFileIsGZipped(outputVCFFileIsGZipped),
	myRefSeq(refSeqFileName.c_str(), false)
{	
	initLogGLMatrix();
}

template <class PILEUP_TYPE, class FUNC_CLASS>
PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::PileupWithGenomeReference(int window, 
																			  const std::string& refSeqFileName,
																			  bool addDelAsBase, 
																			  bool inputVCFFileIsGZipped,
							  												  bool outputVCFFileIsGZipped, 
																			  const FUNC_CLASS& fp)
:	Pileup<PILEUP_TYPE>(window),
	myAddDelAsBase(addDelAsBase),
	inputVCFFileIsGZipped(inputVCFFileIsGZipped),
	outputVCFFileIsGZipped(outputVCFFileIsGZipped),
	myRefSeq(refSeqFileName.c_str()),
	myOutputVCFFile(NULL)
{	
	initLogGLMatrix();
}

template <class PILEUP_TYPE, class FUNC_CLASS>
int PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processFile(const std::string& bamFileName,
																	const std::string& outputVCFFileName, 
                                                					uint16_t excludeFlag,
                                                 					uint16_t includeFlag)
{
	myOutputVCFFile = new InputFile(outputVCFFileName.c_str(), "w", outputVCFFileIsGZipped?InputFile::GZIP:InputFile::DEFAULT);
	if ( !myOutputVCFFile->isOpen() ) {
	  std::cerr << "FATAL error : Cannot open file " << outputVCFFileName << std::endl;
	  abort();
	}
	std::string tempStr("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
	tempStr.append(outputVCFFileName.c_str());
	tempStr.append("\n");
	myOutputVCFFile->ifwrite(tempStr.c_str(), tempStr.length());
	
	Pileup<PILEUP_TYPE>::processFile(bamFileName);
		
	myOutputVCFFile->ifclose();
   
    return(0);  
}

template <class PILEUP_TYPE, class FUNC_CLASS>
int PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processFile(const std::string& bamFileName, 
																	const std::string& bamIndexFileName, 
																    const std::string& inputVCFFileName, 
                                                					const std::string& outputVCFFileName, 
                                                					uint16_t excludeFlag,
                                                 					uint16_t includeFlag)
{
	myOutputVCFFile = new InputFile(outputVCFFileName.c_str(), "w", outputVCFFileIsGZipped?InputFile::GZIP:InputFile::DEFAULT);
	if ( !myOutputVCFFile->isOpen() ) {
	  std::cerr << "FATAL error : Cannot open file " << outputVCFFileName << std::endl;
	  abort();
	}
	std::string tempStr("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
	tempStr.append(outputVCFFileName.c_str());
	tempStr.append("\n");
	myOutputVCFFile->ifwrite(tempStr.c_str(), tempStr.length());
		    	
	//read through input VCF file to collect regions
	InputFile vcfIn(inputVCFFileName.c_str(), "r", inputVCFFileIsGZipped?InputFile::GZIP:InputFile::DEFAULT);
	if ( !vcfIn.isOpen() ) {
	  std::cerr << "FATAL error : Cannot open file " << inputVCFFileName << std::endl;
	  abort();
	}

	std::string field;
	std::string chromosome = "0";
	uint32_t position = 0;
	char c = ' ';
	vector<Region> regions;
	Region currentRegion;
	currentRegion.chrom = "0";
	currentRegion.start = 0;
	currentRegion.end = 0;
	currentRegion.positions = new vector<uint32_t>();
	currentRegion.currentPosition = 0;
	
	while(!vcfIn.ifeof())
	{
		if((c=vcfIn.ifgetc())!='#')
		{
			field.append(&c, 1);
			
			//read chromosome
			while((c=vcfIn.ifgetc())!='\t')
			{
				field.append(&c, 1);
			}
			
			//std::cout << "chrom:" << field << " ";
			chromosome = field;
			field.clear();
			
			//read position
			while((c=vcfIn.ifgetc())!='\t')
			{
				field.append(&c, 1);
			}
			
			//std::cout << " pos:" << field << "\n";
			position = (uint32_t) atoi(field.c_str());		
			field.clear();
			
			//std::cout << "X chrom:" << chromosome << " pos:" << position << "\n";

			//decide to include region or not
			if(chromosome!=currentRegion.chrom)
			{
				//add current region
				Region newRegion;
				newRegion.chrom = currentRegion.chrom;
				newRegion.start = currentRegion.start;
				newRegion.end = currentRegion.end;
				newRegion.positions = currentRegion.positions;
				newRegion.currentPosition = currentRegion.currentPosition;
				regions.push_back(newRegion);
				
				//create new current region
				currentRegion.chrom=chromosome;
				currentRegion.start=position;
				currentRegion.end=position;
				currentRegion.positions = new vector<uint32_t>();
				currentRegion.positions->push_back(position);
				currentRegion.currentPosition = 0;		
			}
			else
			{
				//extend region
				//std::cout << position << " minus " << currentRegion.end << "\n"; 
				if(position-currentRegion.end<1000)
				{
					currentRegion.end = position;
					currentRegion.positions->push_back(position);
				}
				else
				{
					//add current region
					Region newRegion;
					newRegion.chrom = currentRegion.chrom;
					newRegion.start = currentRegion.start;
					newRegion.end = currentRegion.end;
					newRegion.positions = currentRegion.positions;
					regions.push_back(newRegion);
					
					//create new current region
					currentRegion.chrom=chromosome;
					currentRegion.start=position;
					currentRegion.end=position;
					currentRegion.positions = new vector<uint32_t>();
					currentRegion.positions->push_back(position);
					
				}	
			}		
		}
		
		//read rest of line
		while((c=vcfIn.ifgetc())!='\n')
		{
		}
	}

	// add the last region
	if(currentRegion.end-currentRegion.start > 0) {
	  Region newRegion;
	  newRegion.chrom = currentRegion.chrom;
	  newRegion.start = currentRegion.start;
	  newRegion.end = currentRegion.end;
	  newRegion.positions = currentRegion.positions;
	  regions.push_back(newRegion);	
	}
	
	/*
	//std::cout << "size " << regions.size()<<"\n";
	for (uint i=1; i<regions.size(); ++i)
	//for (uint i=86770; i<86771; ++i)
	{
		std::cout << "chr" << regions.at(i).chrom.c_str() << ":" <<  regions.at(i).start << "-" << regions.at(i).end;
		std::cout << " positions: " ;
		
		for (uint k=0; k<regions.at(i).positions->size(); ++k)	
		{
			std::cout << regions.at(i).positions->at(k) << " ";
		}
		
		std::cout  << "\n";
	}
	*/
    SamFile samIn;
    SamFileHeader header;
    SamRecord record;
    
    if(!samIn.OpenForRead(bamFileName.c_str()))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }
    
    if(!samIn.ReadHeader(header))
    {
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        return(samIn.GetStatus());
    }

  	// read index file
  	if (!samIn.ReadBamIndex(bamIndexFileName.c_str()))
  	{
  		//std::cout << "samin status" << samIn.GetStatus() <<"\n";
		return(samIn.GetStatus());  	
	}
	
    // The file needs to be sorted by coordinate.
    samIn.setSortedValidation(SamFile::COORDINATE);

	// Iterate over selected regions
	for (uint i=1; i<regions.size(); ++i)
	//for (uint i=86770; i<86771; ++i)
	{
		int lastReadAlignmentStart = 0;
		Region currentRegion = regions.at(i);
		if(!samIn.SetReadSection(header.getReferenceID(currentRegion.chrom.c_str()), currentRegion.start-1, currentRegion.end))
		{
			std::cout << "chrom:" << regions.at(i).chrom << ":" <<  regions.at(i).start << "-" << regions.at(i).end << "\n";
		}
		else
		{
			//std::cout << "reading chrom:" << regions.at(i).chrom << ":" <<  regions.at(i).start << "-" << regions.at(i).end << "\n";
			
			//int count = 0;
			 // Iterate over all records
		    while (samIn.ReadRecord(header, record))
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
		        
		        if(record.get0BasedPosition()>=lastReadAlignmentStart)
		        {
		        	lastReadAlignmentStart = record.get0BasedPosition();
		        	//std::cout << "lastReadAlignmentStart " << lastReadAlignmentStart << "\n";
		        	processAlignment(record, &currentRegion);
		        }
		        		        	
	        	//++count;
			}
			
			//std::cout << " count:" << count << "\n";			
	    }
	}	
   
    Pileup<PILEUP_TYPE>::flushPileup();

    int returnValue = 0;
    if(samIn.GetStatus() != SamStatus::NO_MORE_RECS)
    {
        // Failed to read a record.
        fprintf(stderr, "%s\n", samIn.GetStatusMessage());
        returnValue = samIn.GetStatus();
    }
 	
	myOutputVCFFile->ifclose();
   
    return(returnValue);  
}

template <class PILEUP_TYPE, class FUNC_CLASS>
void PileupWithGenomeReference<PILEUP_TYPE, FUNC_CLASS>::processAlignment(SamRecord& record, Region* region)
{
    int refPosition = record.get0BasedPosition();
    int refID = record.getReferenceID();

    // Flush any elements from the pileup that are prior to this record
    // since the file is sorted, we are done with those positions.
    //std::cout << "        flushing before " << refID << ":"<< refPosition << "\n";
    Pileup<PILEUP_TYPE>::flushPileup(refID, refPosition);
    
    // Loop through for each reference position covered by the record.
    // It is up to the PILEUP_TYPE to handle insertions/deletions, etc
    // that are related with the given reference position.
    //for(; refPosition <= record.get0BasedAlignmentEnd(); ++refPosition)
    //{
    //    Pileup<PILEUP_TYPE>::addAlignmentPosition(refPosition, record);
    //}
	
	//search for first location in region.positions that is >= the start position of the record
	while (region->positions->at(region->currentPosition)-1< (uint32_t)refPosition)
	{
		++(region->currentPosition);
	}
	//std::cout << "        currentPosition " << region->currentPosition << "\n";
	
	for (uint k=region->currentPosition; k<region->positions->size()&& region->positions->at(k)-1<=(uint32_t)record.get0BasedAlignmentEnd(); ++k)	
	{
		//std::cout << "            adding Position " << region->positions->at(k) << "\n";
		Pileup<PILEUP_TYPE>::addAlignmentPosition(region->positions->at(k)-1, record);
	}
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
