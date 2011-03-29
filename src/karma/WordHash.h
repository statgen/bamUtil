/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#if !defined(WORDHASH_H)
#define WORDHASH_H

#include "GenomeSequence.h"
#include "MemoryMapArray.h"
// #include "WordIndex.h"
#include <map>
#include <vector>

using std::map;
using std::vector;
using std::pair;

#if !defined(HAVE_WORDINTEGER_T)
#define HAVE_WORDINTEGER_T
typedef uint32_t wordInteger_t;
typedef uint32_t wordIndex_t;
#define INVALID_WORDINDEX UINT32_MAX
#endif

//
// There are a pair of hashes - I call them the left and right hash.
// The hash keys are formed by using a high repeat word combined with
// the word following (the right hash), and the same high repeat word
// combined with the word preceding (the left hash).
//
// Each hash key is used as a primary key into a hash table of total
// length 1.33 times the number of unique keys rounded up to the nearest
// prime number.
//
// The hash table entries are composed of a 128-bit structure with: the
// 64 bit hash key, and a 32 bit word position table index, and lastly,
// a 32 bit word count.
//
// In collisions, a secondary prime modulus (1013) is used to skip to
// another table entry until the correct key is found.
//
// Once found, the hash lookup returns the start of the corresponding
// location in the word positions table, and the count, so that generic
// higher level code can evaluate the positions.
//

//
// implement a fast lookup from a pair of wordInteger_t values to
// the word positions which they map to.
//
// Typically, one of the two values will be a wordInteger that
// exceeded the max cutoff count for WordIndex (default 5000).
//
// This has is used to get a fast secondary lookup of match
// candidate positions.
//
// The class is derived from the MemoryMapArray object, and is
// simply a fixed length vector of uint32_t values.  It is
// two tables in succession - the first is the small hash, each
// hash entry consisting of a 3-tuple {key, count, wordPositionStartIndex}.
//
// I expect karma to use a pair of these, one for {HIGH_REPEAT, XXXX} and
// {XXXX, HIGH_REPEAT} combinations, the choice depending on which of the
// first several words of the read are repeats.
//
// The goal here is to improve both sensitivity and speed of karma, which
// tends to bog down in the genome areas with high repeats.
//
//
//    hash the key (~10K values) to a tuple {key, count, wordPositionStartIndex},
//    where wordPositionStartIndex is a set of 'count' values in the wordPositions table.
//

#define UMWH_COOKIE 0x38c59d8f  // unique cookie id
#define UMWH_VERSION 20090713U  // YYYYMMDD of last change to file layout

class WordHashMmapHeader : public MemoryMapArrayHeader
{
 public:
    size_t getHeaderSize(int i)
    {
        return sizeof(*this);
    };

    uint32_t    hashCount;      // how many elements in the hash
    uint32_t    hashSize;       // size (==smallest prime that is >=hashCount)
    uint32_t    hashModulus;    // unsure if necessary
    uint32_t    hashModulus2;   // secondary hash modulus
    uint32_t    wordTableSizeStart;
    uint32_t    wordTableSize;
};

std::ostream &operator << (std::ostream &stream, WordHashMmapHeader &h);

typedef MemoryMapArray<
uint32_t,
    uint32_t,
    UMWH_COOKIE,
    UMWH_VERSION,
    mmapUint32Access,
    mmapUint32Set,
    mmapUint32elementCount2Bytes,
    WordHashMmapHeader
    > WordHashArray;


class WordHash : public WordHashArray
{
 private:
    typedef uint64_t    WordHashKey;

    WordHashKey    hash(wordInteger_t prefix, wordInteger_t suffix)
    {
        return (WordHashKey) prefix << (sizeof(suffix) << 3) | suffix;
    }

    typedef  map<WordHashKey, vector<genomeIndex_t> >    GenomeLocations;
    GenomeLocations genomeLocations;

#define UNUSED_ENTRY UINT64_MAX
    struct hashEntry
    {
        WordHashKey key;
        uint32_t    count;
        uint32_t    wordPositions;
    };
    genomeIndex_t   *wordPositions;
    hashEntry       *wordHash;

 public:
    WordHash();
    ~WordHash();
    void close();

    // to create, first call addGenomeLocation() a whole bunch of times,
    // then call create().
    void addGenomeLocation(wordInteger_t prefix, wordInteger_t suffix, genomeIndex_t index);
    bool create(const char *outputFile);

    // to use, first call open, then call findGenomeLocations()
    // a whole bunch of times.
    int open(const char *inputFile);
    int findGenomeLocations(wordInteger_t prefix, wordInteger_t suffix, genomeIndex_t *&values);
};

inline void WordHash::addGenomeLocation(wordInteger_t prefix, wordInteger_t suffix, genomeIndex_t index)
{
    WordHashKey h = hash(prefix, suffix);

    GenomeLocations::iterator it = genomeLocations.find(h);
    if (it==genomeLocations.end())
    {
        pair<WordHashKey, vector<genomeIndex_t> > p;
        p.first = h;
        genomeLocations.insert(p);
    }
    genomeLocations[h].push_back(index);
}

#endif
