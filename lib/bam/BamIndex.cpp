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

#include "BamIndex.h"
#include <iomanip>

Chunk SortedChunkList::pop()
{
    Chunk newChunk = chunkList.begin()->second;
    chunkList.erase(chunkList.begin());
    return(newChunk);
}


bool SortedChunkList::insert(const Chunk& chunkToInsert)
{
    std::pair<std::map<uint64_t, Chunk>::iterator, bool> insertRes;
    // Insert the passed in chunk.
    insertRes = 
        chunkList.insert(std::pair<uint64_t, Chunk>(chunkToInsert.chunk_beg,
                                                    chunkToInsert));

    if(!insertRes.second)
    {
        // Failed to insert the chunk.
        std::cerr << "Failed to insert into the SortedChunkList.\n";
        std::cerr << "\tpreviously found chunk:\tbeg = " << std::hex
                  << insertRes.first->second.chunk_beg 
                  << "\tend = "
                  << insertRes.first->second.chunk_end
                  << "\nnew chunk:\tbeg = " << std::hex
                  << chunkToInsert.chunk_beg 
                  << "\tend = "
                  << chunkToInsert.chunk_end
                  << std::endl;
    }
    // return the result that comes from insertRes.
    return(insertRes.second);
}

void SortedChunkList::clear()
{
    chunkList.clear();
}

bool SortedChunkList::empty()
{
    return(chunkList.empty());
}


// Merge overlapping chunks found in this list.
bool SortedChunkList::mergeOverlapping()
{
    // Start at the beginning of the list and iterate through.
    std::map<uint64_t, Chunk>::iterator currentPos = chunkList.begin();
    std::map<uint64_t, Chunk>::iterator nextPos = chunkList.begin();
    if(nextPos != chunkList.end())
    {
        ++nextPos;
    }
    
    // Loop until the end is reached.
    while(nextPos != chunkList.end())
    {
        // If the next chunk is completely contained within the current 
        // chunk (its end is less than the current chunk's end), 
        // delete it since its position is already covered.
        if(nextPos->second.chunk_end < currentPos->second.chunk_end)
        {
            chunkList.erase(nextPos);
            nextPos = currentPos;
            ++nextPos;
            continue;
        }
        
        // If the next chunk's start position's BGZF block is less than or
        // equal to the BGZF block of the current chunk's end position,
        // combine the two chunks into the current chunk.
        if((nextPos->second.chunk_beg >> 16) <= 
           (currentPos->second.chunk_end >> 16))
        {
            currentPos->second.chunk_end = nextPos->second.chunk_end;
            // nextPos has now been included in the current pos, so
            // remove it.
            chunkList.erase(nextPos);
            nextPos = currentPos;
            ++nextPos;
            continue;
        }
        else
        {
            // Nothing to combine.  So try combining at the next 
            currentPos = nextPos;
            ++nextPos;
        }
    }
    return(true);
}


BamIndex::BamIndex()
{
    resetIndex();
}



BamIndex::~BamIndex()
{
}


// Reset the member data for a new index file.
void BamIndex::resetIndex()
{
    n_ref = 0;
    maxOverallOffset = 0;
    
    myUnMappedNumReads = -1;

    // Clear the references.
    myRefs.clear();
}


// Read & parse the specified index file.
SamStatus::Status BamIndex::readIndex(const char* filename)
{
    // Reset the index from anything that may previously be set.
    resetIndex();

    IFILE indexFile = ifopen(filename, "rb");

    // Failed to open the index file.
    if(indexFile == NULL)
    {
        return(SamStatus::FAIL_IO);
    }

    // generate the bam index structure.

    // Read the magic string.
    char magic[4];
    if(ifread(indexFile, magic, 4) != 4)
    {
        // Failed to read the magic
        return(SamStatus::FAIL_IO);
    }

    // If this is not an index file, set num references to 0. 
    if (magic[0] != 'B' || magic[1] != 'A' || magic[2] != 'I' || magic[3] != 1)
    {
        // Not a BAM Index file.
        return(SamStatus::FAIL_PARSE);
    }

    // It is a bam index file.
    // Read the number of reference sequences.
    if(ifread(indexFile, &n_ref, 4) != 4)
    {
        // Failed to read.
        return(SamStatus::FAIL_IO);
    }

    // Size the references.
    myRefs.resize(n_ref);

    for(int refIndex = 0; refIndex < n_ref; refIndex++)
    {
        // Read each reference.
        Reference* ref = &(myRefs[refIndex]);
        
        // Resize the bins so they can be indexed by bin number.
        ref->bins.resize(MAX_NUM_BINS + 1);
        
        // Read the number of bins.
        if(ifread(indexFile, &(ref->n_bin), 4) != 4)
        {
            // Failed to read the number of bins.
            // Return failure.
            return(SamStatus::FAIL_PARSE);
        }

        // If there are no bins, then there are no
        // mapped/unmapped reads.
        if(ref->n_bin == 0)
        {
            ref->n_mapped = 0;
            ref->n_unmapped = 0;
        }

        // Read each bin.
        for(int binIndex = 0; binIndex < ref->n_bin; binIndex++)
        {
            uint32_t binNumber;

            // Read in the bin number.
            if(ifread(indexFile, &(binNumber), 4) != 4)
            {
                // Failed to read the bin number.
                // Return failure.
                return(SamStatus::FAIL_IO);
            }

            // Add the bin to the reference and get the
            // pointer back so the values can be set in it.
            Bin* binPtr = &(ref->bins[binNumber]);
            binPtr->bin = binNumber;
         
            // Read in the number of chunks.
            if(ifread(indexFile, &(binPtr->n_chunk), 4) != 4)
            {
                // Failed to read number of chunks.
                // Return failure.
                return(SamStatus::FAIL_IO);
            }

            // Read in the chunks.
            // Allocate space for the chunks.
            uint32_t sizeOfChunkList = binPtr->n_chunk * sizeof(Chunk);
            binPtr->chunks = (Chunk*)malloc(sizeOfChunkList);
            if(ifread(indexFile, binPtr->chunks, sizeOfChunkList) != sizeOfChunkList)
            {
                // Failed to read the chunks.
                // Return failure.
                return(SamStatus::FAIL_IO);
            }

            // Determine the min/max for this bin if it is not the max bin.
            if(binPtr->bin != MAX_NUM_BINS)
            {
                for(int i = 0; i < binPtr->n_chunk; i++)
                {
                    if(binPtr->chunks[i].chunk_beg < ref->minChunkOffset)
                    {
                        ref->minChunkOffset = binPtr->chunks[i].chunk_beg;
                    }
                    if(binPtr->chunks[i].chunk_end > ref->maxChunkOffset)
                    {
                        ref->maxChunkOffset = binPtr->chunks[i].chunk_end;
                    }
                    if(binPtr->chunks[i].chunk_end > maxOverallOffset)
                    {
                        maxOverallOffset = binPtr->chunks[i].chunk_end;
                    }
                }
            }
            else
            {
                // Mapped/unmapped are the last chunk of the
                // MAX BIN
                ref->n_mapped = binPtr->chunks[binPtr->n_chunk - 1].chunk_beg;
                ref->n_unmapped = binPtr->chunks[binPtr->n_chunk - 1].chunk_end;
            }
        }

        // Read the number of intervals.
        if(ifread(indexFile, &(ref->n_intv), 4) != 4)
        {
            // Failed to read, set to 0.
            ref->n_intv = 0;
            // Return failure.
            return(SamStatus::FAIL_IO);
        }

        // Allocate space for the intervals and read them.
        uint32_t linearIndexSize = ref->n_intv * sizeof(uint64_t);
        ref->ioffsets = (uint64_t*)malloc(linearIndexSize);
        if(ifread(indexFile, ref->ioffsets, linearIndexSize) != linearIndexSize)
        {
            // Failed to read the linear index.
            // Return failure.
            return(SamStatus::FAIL_IO);
        }
    }

    int32_t numUnmapped = 0;
    if(ifread(indexFile, &numUnmapped, sizeof(int32_t)) == sizeof(int32_t))
    {
        myUnMappedNumReads = numUnmapped;
    }

    // Successfully read teh bam index file.
    return(SamStatus::SUCCESS);
}


// Get the chunks for the specified reference id and start/end 0-based
// coordinates.
bool BamIndex::getChunksForRegion(int32_t refID, int32_t start, int32_t end, 
                                  SortedChunkList& chunkList)
{
    chunkList.clear();

    // If start is >= to end, there will be no sections, return no
    // regions.
    if((start >= end) && (end != -1))
    {
        std::cerr << "Warning, requesting region where start <= end, so "
                  << "no values will be returned.\n";
        return(false);
    }

    // Handle REF_ID_UNMAPPED.  This uses a default chunk which covers
    // from the max offset to the end of the file.
    if(refID == REF_ID_UNMAPPED)
    {
        Chunk refChunk;
        // The start of the unmapped region is the max offset found
        // in the index file.
        refChunk.chunk_beg = getMaxOffset();
        // The end of the unmapped region is the end of the file, so
        // set chunk end to the max value.
        refChunk.chunk_end = Chunk::MAX_CHUNK_VALUE;
        return(chunkList.insert(refChunk));
    }

    if((refID < 0) || (refID >= n_ref))
    {
        // The specified refID is out of range, return false.
        return(false);
    }

    const Reference* ref = &(myRefs[refID]);

    // Handle where start/end are defaults.    
    if(start == -1)
    {
        if(end == -1)
        {
            // This is whole chromosome, so take a shortcut.
            if(ref->maxChunkOffset == 0)
            {
                // No chunks for this region, but this is not an error.
                return(true);
            }
            Chunk refChunk;
            refChunk.chunk_beg = ref->minChunkOffset;
            refChunk.chunk_end = ref->maxChunkOffset;
            return(chunkList.insert(refChunk));
        }
        else
        {
            start = 0;
        }
    }
    if(end == -1)
    {
        // MAX_POSITION is inclusive, but end is exclusive, so add 1.
        end = MAX_POSITION + 1;
    }

    // Determine the minimum offset for the given start position.  This
    // is done by using the linear index for the specified start position.
    uint64_t minOffset = getMinOffsetFromLinearIndex(refID, start);

    uint16_t binInRangeList[MAX_NUM_BINS + 1];
    
    int numBins = getBinsForRegion(start, end, binInRangeList);

    // loop through the bins in the range and get the chunks.
    for(int i = 0; i < numBins; ++i)
    {
        int binNum = binInRangeList[i];
        const Bin* bin = &(ref->bins[binNum]);

        // Add each chunk in the bin to the map.
        for(int chunkIndex = 0; chunkIndex < bin->n_chunk; chunkIndex++)
        {
            // If the end of the chunk is less than the minimum offset
            // for the 16K block that starts our region, then no
            // records in this chunk will cross our region, so do
            // not add it to the chunks we need to use.
            if(bin->chunks[chunkIndex].chunk_end < minOffset)
            {
                continue;
            }
            // Add the chunk to the map.
            if(!chunkList.insert(bin->chunks[chunkIndex]))
            {
                // Failed to add to the map, return false.
                return(false);
            }
        }
    }

    // Now that all chunks have been added to the list,
    // handle overlapping chunks.
    return(chunkList.mergeOverlapping());
}


// Get the max offset.
uint64_t BamIndex::getMaxOffset() const
{
    return(maxOverallOffset);
}

// Get the min & max file offsets for the reference ID.
bool BamIndex::getReferenceMinMax(int32_t refID,
                                  uint64_t& minOffset,
                                  uint64_t& maxOffset) const
{
    if((refID < 0) || (refID >= (int32_t)myRefs.size()))
    {
        // Reference ID is out of range for this index file.
        return(false);
    }

    // Get this reference.
    minOffset = myRefs[refID].minChunkOffset;
    maxOffset = myRefs[refID].maxChunkOffset;
    return(true);
}

    
// Get the number of references in this index.
int32_t BamIndex::getNumRefs() const
{
    // Return the number of references.
    return(myRefs.size());
}


// Get the number of mapped reads for this reference id.
int32_t BamIndex::getNumMappedReads(int32_t refID)
{
    // If it is the reference id of unmapped reads, return
    // that there are no mapped reads.
    if(refID == REF_ID_UNMAPPED)
   {
       // These are by definition all unmapped reads.
       return(0);
   }

    if((refID < 0) || (refID >= (int32_t)myRefs.size()))
    {
        // Reference ID is out of range for this index file.
        return(-1);
    }

    // Get this reference.
    return(myRefs[refID].n_mapped);
}


// Get the number of unmapped reads for this reference id.
int32_t BamIndex::getNumUnMappedReads(int32_t refID)
{
    // If it is the reference id of unmapped reads, return
    // that value.
    if(refID == REF_ID_UNMAPPED)
    {
        return(myUnMappedNumReads);
    }

    if((refID < 0) || (refID >= (int32_t)myRefs.size()))
    {
        // Reference ID is out of range for this index file.
        return(-1);
    }

    // Get this reference.
    return(myRefs[refID].n_unmapped);
}


// Print the bam index.
void BamIndex::printIndex(int32_t refID, bool summary)
{
    std::cout << "BAM Index: " << std::endl;
    std::cout << "# Reference Sequences: " << n_ref << std::endl;

    unsigned int startRef = 0;
    unsigned int endRef = myRefs.size() - 1;
    std::vector<Reference> refsToProcess;
    if(refID != -1)
    {
        // Set start and end ref to the specified reference id.
        startRef = refID;
        endRef = refID;
    }

    // Print out the information for each bin.
    for(unsigned int i = startRef; i <= endRef; ++i)
    {
        std::cout << std::dec 
                  << "\tReference ID: " << std::setw(4) << i
                  << ";  #Bins: "<< std::setw(6) << myRefs[i].n_bin 
                  << ";  #Linear Index Entries: " 
                  << std::setw(6) << myRefs[i].n_intv
                  << ";  Min Chunk Offset: " 
                  << std::setw(18) << std::hex << std::showbase << myRefs[i].minChunkOffset
                  << ";  Max Chunk Offset: "
                  << std::setw(18) << myRefs[i].maxChunkOffset
                  << std::dec;
        // Print the mapped/unmapped if set.
        if(myRefs[i].n_mapped != Reference::UNKNOWN_MAP_INFO)
        {            
            std::cout << ";  " << myRefs[i].n_mapped << " Mapped Reads";
        }
        if(myRefs[i].n_mapped != Reference::UNKNOWN_MAP_INFO)
        {            
            std::cout << ";  " << myRefs[i].n_unmapped << " Unmapped Reads";
        }
        std::cout << std::endl;
        
        // Only print more details if not summary.
        if(!summary)
        {
            std::vector<Bin>::iterator binIter;
            for(binIter = myRefs[i].bins.begin(); 
                binIter != myRefs[i].bins.end();
                ++binIter)
            {
                Bin* binPtr = &(*binIter);
                if(binPtr->bin == Bin::NOT_USED_BIN)
                {
                    // This bin is not used, continue.
                    continue;
                }
                // Print the bin info.
                std::cout << "\t\t\tBin Name: " << binPtr->bin << std::endl;
                std::cout << "\t\t\t# Chunks: " << binPtr->n_chunk << std::endl;
                std::cout << std::hex << std::showbase;

                for(int chunkIndex = 0; chunkIndex < binPtr->n_chunk;
                    ++chunkIndex)
                {
                    // If this is the last chunk of the MAX_NUM_BINS - it
                    // contains a mapped/unmapped count rather than the regular
                    // chunk addresses.
                    if((binPtr->bin != MAX_NUM_BINS) ||
                       (chunkIndex != (binPtr->n_chunk - 1)))
                    {
                        std::cout << "\t\t\t\tchunk_beg: "
                                  << binPtr->chunks[chunkIndex].chunk_beg 
                                  << std::endl;
                        std::cout << "\t\t\t\tchunk_end: "
                                  << binPtr->chunks[chunkIndex].chunk_end
                                  << std::endl;
                    }
                }
            }
            std::cout << std::dec;
            
            // Print the linear index.
            for(int linearIndex = 0; linearIndex < myRefs[i].n_intv;
                ++linearIndex)
            {
                if(myRefs[i].ioffsets[linearIndex] != 0)
                {
                    std::cout << "\t\t\tLinearIndex["
                              << std::dec << linearIndex << "] Offset: " 
                              << std::hex << myRefs[i].ioffsets[linearIndex]
                              << std::endl;
                }
            }
        }
    }
}


// Returns the minimum offset of records that cross the 16K block that
// contains the specified position for the given reference id.
uint64_t BamIndex::getMinOffsetFromLinearIndex(int32_t refID, uint32_t position) const
{
    int32_t linearIndex = position >> LINEAR_INDEX_SHIFT;
    uint64_t minOffset = 0;
    // Check to see if the position is out of range of the linear index.
    int32_t linearOffsetSize = myRefs[refID].n_intv;

    // If there are no entries in the linear index, return 0.
    // Or if the linear index is not large enough to include
    // the start block, then there can be no records that cross
    // our region, so return 0.
    if((linearOffsetSize == 0) || (linearIndex >= linearOffsetSize))

    {
        return(0);
    }
    else
    {
        // The linear index is specified for this block, so return that
        // value.
        minOffset = myRefs[refID].ioffsets[linearIndex];

        // If the offset is 0, go to the previous block that has an offset.
        // This is due to a couple of bugs in older sam tools indexes.
        // 1) they add one to the index location (so when reading those, you
        // may be starting earlier than necessary)
        // 2) (the bigger issue) They did not include bins 4681-37449 in
        // the linear index.
        while((minOffset == 0) && (--linearIndex >= 0))
        {
            minOffset = myRefs[refID].ioffsets[linearIndex];            
        }
    }
    return(minOffset);
}


// Returns the minimum offset of records that cross the 16K block that
// contains the specified position for the given reference id.
// The basic logic is from samtools reg2bins and the samtools format specification pdf.
int BamIndex::getBinsForRegion(uint32_t start, uint32_t end, uint16_t binList[MAX_NUM_BINS])
{
	uint32_t binListIndex = 0, binNum;
	--end;

        // Check if beg/end go too high, set to max position.
        if(start > MAX_POSITION)
        {
            start = MAX_POSITION;
        }
        if(end > MAX_POSITION)
        {
            end = MAX_POSITION;
        }

	binList[binListIndex++] = 0;
	for (binNum =    1 + (start>>26); binNum <=    1 + (end>>26); ++binNum) 
            binList[binListIndex++] = binNum;
	for (binNum =    9 + (start>>23); binNum <=    9 + (end>>23); ++binNum) 
            binList[binListIndex++] = binNum;
	for (binNum =   73 + (start>>20); binNum <=   73 + (end>>20); ++binNum)
            binList[binListIndex++] = binNum;
	for (binNum =  585 + (start>>17); binNum <=  585 + (end>>17); ++binNum)
            binList[binListIndex++] = binNum;
	for (binNum = 4681 + (start>>14); binNum <= 4681 + (end>>14); ++binNum)
            binList[binListIndex++] = binNum;

        // binListIndex contains the number of items added to the list.
	return binListIndex;
}


