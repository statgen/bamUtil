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

#if !defined(_CIGAR_H)
#define _CIGAR_H

#include <string.h> // for inline use of strcat, etc
#include <limits.h> // for INT_MAX
#include <stdint.h> // for uint32_t and friends


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>

#include "Generic.h"
#include "StringBasics.h"

//
// Docs from Sam1.pdf:
//
// Clipped alignment. In Smith-Waterman alignment, a sequence may not be aligned from the first residue to the last one.
// Subsequences at the ends may be clipped off. We introduce operation ʻSʼ to describe (softly) clipped alignment. Here is
// an example. Suppose the clipped alignment is:
// REF:  AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG
// READ:        gggGTGTAACC-GACTAGgggg
// where on the read sequence, bases in uppercase are matches and bases in lowercase are clipped off. The CIGAR for
// this alignment is: 3S8M1D6M4S.
//
//
// If the mapping position of the query is not available, RNAME and
// CIGAR are set as “*”
//
// A CIGAR string is comprised of a series of operation lengths plus the operations. The conventional CIGAR format allows
// for three types of operations: M for match or mismatch, I for insertion and D for deletion. The extended CIGAR format
// further allows four more operations, as is shown in the following table, to describe clipping, padding and splicing:
//
// op   Description
// --   -----------
// M    Match or mismatch
// I    Insertion to the reference
// D    Deletion from the reference
// N    Skipped region from the reference
// S    Soft clip on the read (clipped sequence present in <seq>)
// H    Hard clip on the read (clipped sequence NOT present in <seq>)
// P    Padding (silent deletion from the padded reference sequence)
//
//



////////////////////////////////////////////////////////////////////////
//
// This class represents the CIGAR.  It contains methods for converting
// to strings and extracting information from the cigar on how a read
// maps to the reference.
//
// It only contains read only methods.  There are no ways to set
// values.  To set a value, a child class must be used.
//
class Cigar
{
public:
    enum Operation {none, match, mismatch, insert, del, skip, softClip, hardClip, pad};

    // Return true if the specified operation is found in the
    // query sequence, false if not.
    static bool foundInQuery(Operation op)
    {
        switch(op)
        {
            case match:
            case mismatch:
            case insert:
            case softClip:
                return true;
            default:
                return false;
        }
        return false;
    }
    
    // Return true if the specified operation is a clipping operation,
    // false if not.
    static bool isClip(Operation op)
    {
        switch(op)
        {
            case softClip:
            case hardClip:
                return true;
            default:
                return false;
        }
        return false;
    }
    
    ////////////////////////////////////////////////////////////////////////
    //
    // Nested Struct : CigarOperator
    //
    struct CigarOperator
    {

        CigarOperator()
        {
            operation = none;
            count = 0;
        }

        CigarOperator(Operation operation, uint32_t count): operation(operation), count(count) {};

        Operation operation;

        uint32_t count;

        // Get the character associated with this operation.
        char getChar() const
        {
            switch (operation)
            {
                case none:
                    return '?';  // error
                case match:
                case mismatch:
                    return'M';
                case insert:
                    return 'I';
                case del:
                    return'D';
                case skip:
                    return 'N';
                case softClip:
                    return 'S';
                case hardClip:
                    return 'H';
                case pad:
                    return 'P';
            }
            return '?'; // actually it is an error to get here
        }

        //
        // Compare only on the operator.
        //
        // Match and mismatch are considered the same for CIGAR strings.
        //
        bool operator == (const CigarOperator &rhs) const
        {
            if (operation==rhs.operation)
                return true;
            if ((operation == mismatch || operation == match) && (rhs.operation == mismatch || rhs.operation == match))
                return true;
            return false;
        }

        bool operator != (const CigarOperator &rhs) const
        {
            return !((*this) == rhs) ;
        }

    };

    ////////////////////////////////////////////////////////////////////////
    //
    // Cigar  Class
    //
    friend std::ostream &operator << (std::ostream &stream, const Cigar& cigar);

    Cigar()
    {
        clearQueryAndReferenceIndexes();
    }

    //
    // Set the passed in string to the string reprentation of the Cigar
    // operations in this object.
    //
    void getCigarString(String& cigarString) const;
    void getCigarString(std::string& cigarString) const;

    /// obtain a non-run length encoded string of operations
    ///
    /// The returned string is actually also a valid CIGAR string,
    /// but it does not have any digits in it - only the characters
    /// themselves.  In theory this makes it easier to parse some
    /// reads.
    ///
    /// /return s the string to populate
    void getExpandedString(std::string &s) const;

    const CigarOperator & operator [](int i) const
    {
        return cigarOperations[i];
    }

    const CigarOperator & getOperator(int i) const
    {
        return cigarOperations[i];
    }

    bool operator == (Cigar &rhs) const;

    //
    // get the number of cigar operations
    //
    int size()  const
    {
        return cigarOperations.size();
    }

    void Dump() const
    {
        String cigarString;
        getCigarString(cigarString);
        std::cout << cigarString ;
    }

    /// return the length of the read that corresponds to
    /// the current CIGAR string.
    ///
    /// For validation, we should expect that a sequence
    /// read in a SAM file will be the same length as the
    /// value returned by this method.
    ///
    /// Example: 3M2D3M describes a read with three bases
    /// matching the reference, then skips 2 bases, then has
    /// three more bases that match the reference (match/mismatch).
    /// In this case, the read length is expected to be 6.
    ///
    /// Example: 3M2I3M describes a read with 3 match/mismatch
    /// bases, two extra bases, and then 3 more match/mistmatch
    /// bases.  The total in this example is 8 bases.
    ///
    /// /return returns the expected read length
    int getExpectedQueryBaseCount() const;

    /// return the number of bases in the reference that
    /// this read "spans"
    ///
    /// When doing range checking, we occassionally need to know
    /// how many total bases the CIGAR string represents as compared
    /// to the reference.
    ///
    /// Examples: 3M2D3M describes a read that overlays 8 bases in
    /// the reference.  3M2I3M describes a read with 3 bases that
    /// match the reference, two additional bases that aren't in the
    /// reference, and 3 more bases that match the reference, so it
    /// spans 6 bases in the reference.
    ///
    /// /return how many bases in the reference are spanned
    /// by the given CIGAR string
    ///
    int getExpectedReferenceBaseCount() const;

    /// Return the number of clips that are at the beginning of the cigar.
    int getNumBeginClips() const;

    /// Return the number of clips that are at the end of the cigar.
    int getNumEndClips() const;

    // Return the reference offset associated with the specified
    // query index based on this cigar.
    int32_t getRefOffset(int32_t queryIndex);

    // Return the query index associated with the specified
    // reference offset based on this cigar.
    int32_t getQueryIndex(int32_t refOffset);

    // Return the reference position associated with the specified query
    // index based on this cigar and the specified queryStartPos which
    // is the leftmost mapping position of the first matching
    // base in the query.
    int32_t getRefPosition(int32_t queryIndex, int32_t queryStartPos);

    // Return the query index associated with the specified reference position
    // when the query starts at the specified reference position based on
    // this cigar.
    int32_t getQueryIndex(int32_t refPosition, int32_t queryStartPos);

    // Return the number of bases that overlap the reference and the
    // read associated with this cigar that falls within the specified region.
    // start : inclusive 0-based start position (reference position) of the
    //         region to check for overlaps in.
    //         (-1 indicates to start at the beginning of the reference.)
    // end   : exclusive 0-based end position (reference position) of the
    //          region to check for overlaps in.
    //         (-1 indicates to go to the end of the reference.)
    // queryStartPos : 0-based leftmost mapping position of the first matching
    //                 base in the query.
    uint32_t getNumOverlaps(int32_t start, int32_t end, int32_t queryStartPos);

    static const int32_t INDEX_NA;

protected:
    // Clear the query index/reference offset index vectors.
    void clearQueryAndReferenceIndexes();

    // Set the query index/reference offset index vectors.
    void setQueryAndReferenceIndexes();

    // Container for the cigar operations in this cigar.
    std::vector<CigarOperator> cigarOperations;

private:
    // The vector is indexed by query index and contains the reference
    // offset associated with that query index.
    // The vector is reset each time a new cigar operation is added, and
    // is calculated when accessed if it is not already set.
    std::vector<int32_t> queryToRef;

    // The vector is indexed by reference offset and contains the query
    // index associated with that reference offset.
    // The vector is reset each time a new cigar operation is added, and
    // is calculated when accessed if it is not already set.
    std::vector<int32_t> refToQuery;
};

inline std::ostream &operator << (std::ostream &stream, const Cigar::CigarOperator& o)
{
    stream << o.count << o.getChar();
    return stream;
}

inline std::ostream &operator << (std::ostream &stream, const Cigar& cigar)
{
    stream << cigar.cigarOperations;
    return stream;
}

#endif
