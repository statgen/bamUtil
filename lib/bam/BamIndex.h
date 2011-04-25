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

#ifndef __BAM_INDEX_H__
#define __BAM_INDEX_H__

#include <stdint.h>
#include <vector>
#include <map>
#include <stdlib.h>

#include "InputFile.h"
#include "SamStatus.h"

class Chunk
{
public:
    uint64_t chunk_beg; // offset of the start of the chunk
    uint64_t chunk_end; // offset of the end of the chunk
    
    static const uint64_t MAX_CHUNK_VALUE = 0xFFFFFFFFFFFFFFFFULL;

    bool operator< (const Chunk& otherChunk) const
    {
        return(this->chunk_beg < otherChunk.chunk_beg);
    }
};


// This class contains chunks that are sorted by the beginning position.
// This class hides how the chunks are actually stored (map, list ,etc),
// so they can be interchanged.
class SortedChunkList
{
public:
    // Returns the first chunk in the list and  removes it.
    Chunk pop();
    bool insert(const Chunk& chunkToInsert);
    void clear();
    bool empty();
    bool mergeOverlapping();

private:
    std::map<uint64_t, Chunk> chunkList;
};

class BamIndex
{
public:

    BamIndex();
    ~BamIndex();

    /// Reset the member data for a new index file.
    void resetIndex();

    // Read & parse the specified index file.
    /// \param filename the bam index file to be read.
    /// \return the status of the read.
    SamStatus::Status readIndex(const char* filename);

    /// Get the list of chunks associated with this region.
    /// For an entire reference ID, set start and end to -1.
    /// To start at the beginning of the region, set start to 0/-1.
    /// To go to the end of the region, set end to -1.
    bool getChunksForRegion(int32_t refID, int32_t start, int32_t end, 
                            SortedChunkList& chunkList);

    uint64_t getMaxOffset() const;

    /// Get the minimum and maximum file offsets for the specfied reference ID.
    /// \param refID the reference ID to locate in the file.
    /// \param minOffset returns the min file offset for the specified reference
    /// \param maxOffset returns the max file offset for the specified reference
    /// \return whether or not the reference was found in the file
    bool getReferenceMinMax(int32_t refID, 
                            uint64_t& minOffset, 
                            uint64_t& maxOffset) const;
    
    /// Get the number of references in this index.
    /// \return number of references
    int32_t getNumRefs() const;

    /// Get the number of mapped reads for this reference id.  Returns -1 for
    /// out of range refIDs.
    /// \param refID reference ID for which to extract the number of mapped reads.
    /// \return number of mapped reads for the specified reference id.
    int32_t getNumMappedReads(int32_t refID);

    /// Get the number of unmapped reads for this reference id.  Returns -1 for
    /// out of range refIDs.
    /// \param refID reference ID for which to extract the number of unmapped reads.
    /// \return number of unmapped reads for the specified reference id
    int32_t getNumUnMappedReads(int32_t refID);

    /// Print the index information.
    /// \param refID reference ID for which to print info for.  -1 means print for all references.
    /// \param summary whether or not to just print a summary (defaults to false).  The summary just contains summary info for each reference and not every bin/chunk.
    void printIndex(int32_t refID, bool summary = false);

    // Returns the minimum offset of records that cross the 16K block that
    // contains the specified position for the given reference id.
    uint64_t getMinOffsetFromLinearIndex(int32_t refID, uint32_t position) const;

    // Number of reference sequences.
    /// The number used for an unknown number of reads.
    static const int32_t UNKNOWN_NUM_READS = -1;

    /// The number used for the reference id of unmapped reads.
    static const int32_t REF_ID_UNMAPPED = -1;

    /// The number used to indicate that all reference ids should be used.
    static const int32_t REF_ID_ALL = -2;

private:

    const static uint32_t MAX_NUM_BINS = 37450; // per specs, at most 37450 bins.
    // Maximum allowed position (inclusive 512MB - 1)
    const static uint32_t MAX_POSITION = 536870911;

    // Number of bits in 1 linear index - how much to shift a position by
    // to determine which offset into the linear index to look for it.
    const static uint32_t LINEAR_INDEX_SHIFT = 14;

    class Bin
    {
    public:
        Bin(){chunks = NULL; reset();}
        ~Bin() {reset();}
        void reset()
        {
            if(chunks != NULL)
            {
                free(chunks);
                chunks = NULL;
            }
            n_chunk = 0; 
            bin = NOT_USED_BIN;
        }
        uint32_t bin; // The bin id.
        int32_t n_chunk; // The number of chunks.
        Chunk* chunks; // The chunks for this bin.
        static const uint32_t NOT_USED_BIN = 0xFFFFFFFF;
    };

    class Reference
    {
        // Add one to the max since there may now be an extra bin containing
        // the mapped/unmapped counts.
    public:
        static const int32_t UNKNOWN_MAP_INFO = -1;
        Reference(){ioffsets = NULL; reset();}
        ~Reference(){reset();}
        void reset()
        { 
            bins.clear(); 
            if(ioffsets != NULL)
            {
                free(ioffsets);
                ioffsets = NULL;
            }
            n_bin = 0; 
            n_intv = 0;
            minChunkOffset = -1;
            maxChunkOffset = 0;
            n_mapped = UNKNOWN_MAP_INFO;
            n_unmapped = UNKNOWN_MAP_INFO;
        }
        int32_t n_bin; // The number of bins.
        int32_t n_intv; // Number of intervals.
        std::vector<Bin> bins;  // The bins for this reference.
        uint64_t* ioffsets; // Offsets of intervals first alignments
        uint64_t minChunkOffset;
        uint64_t maxChunkOffset;
        int32_t n_mapped; // Number of mapped reads.
        int32_t n_unmapped; // Number of unmapped reads.
    };

    // Add the bins associated with the specified region to the passed in list.
    // start is incluive, end is exclusive.
    static int getBinsForRegion(uint32_t start, uint32_t end, uint16_t binList[MAX_NUM_BINS]);

    // Number of reference sequences.
    int32_t n_ref;

    uint64_t maxOverallOffset;

    int32_t myUnMappedNumReads;

    // The references.
    std::vector<Reference> myRefs;
};


#endif
