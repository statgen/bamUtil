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

#include "Performance.h"
#include "WordHash.h"
#include "Util.h"       // for ClosestPrimeAbove


WordHash::WordHash()
{
    wordPositions = NULL;
    wordHash = NULL;
}

WordHash::~WordHash()
{
    close();
}

void WordHash::close()
{
    WordHashArray::close();
    genomeLocations.clear();
}

int WordHash::findGenomeLocations(wordInteger_t prefix, wordInteger_t suffix, genomeIndex_t *&values)
{
    WordHashKey h = hash(prefix, suffix);

    //
    // genomeLocations will contain elements if we have been building the
    // hash in this run.  Otherwise, we've done a memory map open, and
    // genomeLocations will have zero elements in it.
    //
    if (genomeLocations.size()>0)
    {
        GenomeLocations::iterator it = genomeLocations.find(h);
        if (it==genomeLocations.end())
        {
            return 0;
        }
        values = &((*it).second.front());
        return (*it).second.size();
    }

    //
    // more normal memory map case:
    //

    uint32_t hashLocation = h % header->hashSize;

    while (wordHash[hashLocation].key != h)
    {
        if (wordHash[hashLocation].key == UNUSED_ENTRY) return 0;
        hashLocation += header->hashModulus2;
        while (hashLocation >= header->hashSize) hashLocation -= header->hashSize;
    }

    values = &(wordPositions[wordHash[hashLocation].wordPositions]);
    return wordHash[hashLocation].count;
}

//
// given a populated genomeLocations map, fill in the memory
// mapped array with the hash of 3-tuples {key, count, wordPositionStartIndex}
// and the wordPositions array.
//
bool WordHash::create(const char *outputFile)
{
    //
    //
    if (genomeLocations.size() == 0)
    {
        //
        // This is not a bug.  It happens when using small genomes
        // where the index is sufficient to cover all locations with
        // a small enough repeats.  Our test case phiX is a perfect
        // example (cutoff 5000, wordSize 11).
        //
        std::cout << "Secondary hash is empty - creating them anyway (typical of small genome).\n";
    }

    GenomeLocations::iterator it;

    Timing t;
    std::cout << "genomeLocations.size() = " << genomeLocations.size() << std::endl;
    std::cout << "add 20%% and 5 elements to get: " << (uint64_t)(genomeLocations.size() * 1.333 + 5) << std::endl;
    std::cout << "computing closest prime above: " << std::flush;
    t.start();

    int closestPrimeAbove = ClosestPrimeAbove(genomeLocations.size() * 1.333 + 5);

    std::cout << closestPrimeAbove << std::endl;

    uint64_t    numberWordPositions = 0;

    for (it = genomeLocations.begin() ; it != genomeLocations.end(); it++)
    {

        assert((*it).second.size()!=0);

        numberWordPositions += (*it).second.size();    // add in length of each vector of word positions
    }

    if (WordHashArray::create(outputFile,
                              closestPrimeAbove *(sizeof(hashEntry) / sizeof(uint32_t)) + numberWordPositions)) return true;

    header->hashCount = genomeLocations.size();
    header->hashSize = closestPrimeAbove;
    header->hashModulus = closestPrimeAbove;
    header->hashModulus2 = 1013;        // good enough for now
    header->wordTableSizeStart = closestPrimeAbove;
    header->wordTableSize = numberWordPositions;
    strncpy(header->application, "Karma Word Hash", sizeof(header->application));

    wordHash = (hashEntry *) data;
    wordPositions = (uint32_t *) data + header->hashSize * (sizeof(hashEntry) / sizeof(uint32_t));

    // mark hash as completely free
    for (unsigned int i=0; i< header->hashSize; i++)
    {
        wordHash[i].key = UNUSED_ENTRY;
    }

    uint64_t wordPosition = 0;
    for (it = genomeLocations.begin() ; it != genomeLocations.end(); it++)
    {
        uint32_t    hashLocation = (*it).first % header->hashModulus;

        while (wordHash[hashLocation].key != UNUSED_ENTRY)
        {
            hashLocation += header->hashModulus2;
            while (hashLocation >= header->hashSize) hashLocation -= header->hashSize;
        }

        wordHash[hashLocation].key = (*it).first;

        wordHash[hashLocation].count = (*it).second.size();

        wordHash[hashLocation].wordPositions = wordPosition;

        for (vector<uint32_t>::iterator wp = (*it).second.begin(); wp != (*it).second.end(); wp++)
            wordPositions[wordPosition++] = *wp;

        numberWordPositions += (*it).second.size();    // add in length of each vector of word positions
    }

    t.end();

    std::cout << "Generated word hash in " << t.interval()/60.0 << " minutes." << std::endl;

    return false;
}

int WordHash::open(const char *inputFile)
{
    if (WordHashArray::open(inputFile)) return true;
    wordHash = (hashEntry *) data;
    // we get away with the division trick here because we are only putting
    // even multiples of uint32_t objects into the hashEntry struct.
    wordPositions = (uint32_t *) data + header->hashSize * (sizeof(hashEntry) / sizeof(uint32_t));
    return false;
}

std::ostream &operator << (std::ostream &stream, WordHashMmapHeader &h)
{
    stream << "hashCount: " << h.hashCount << std::endl;
    stream << "hashSize: " << h.hashSize << std::endl;
    stream << "hashModulus: " << h.hashModulus << std::endl;
    stream << "hashModulus2: " << h.hashModulus2 << std::endl;
    stream << "wordTableSize: " << h.wordTableSize<< std::endl;
    stream << "wordTableSizeStart: " << h.wordTableSizeStart << std::endl;
    return stream;
}




#if defined(TEST)
//
// compile test with:
//
// g++ -g -o testWordHash -DTEST WordHash.cpp Util.o -I../libcsg ../libcsg/libcsg.a -lz
//
int main(int argc, char **argv)
{
    WordHash wh;
    int count;
    genomeIndex_t *positions;

    wh.addGenomeLocation(1,2,3);

    wh.addGenomeLocation(3,4,5);
    wh.addGenomeLocation(3,4,10);
    wh.addGenomeLocation(3,4,100);

    wh.addGenomeLocation(7,1,100);
    wh.addGenomeLocation(7,1,101);
    wh.addGenomeLocation(7,1,102);

    wh.addGenomeLocation(555933597, 927287928, 687443958);
    wh.addGenomeLocation(555933597, 1073741816, 2691223814);
    wh.addGenomeLocation(555933597, 927287928, 2691223814);
    wh.addGenomeLocation(555933597, 927287928, 2874712865);

    count = wh.findGenomeLocations(3,4,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    count = wh.findGenomeLocations(1,2,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    count = wh.findGenomeLocations(7,1,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    count = wh.findGenomeLocations(555933597,927287928,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    wh.create("tooty.hash");

    wh.close();

    wh.open("tooty.hash");

    count = wh.findGenomeLocations(3,4,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    count = wh.findGenomeLocations(1,2,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    count = wh.findGenomeLocations(7,1,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";

    count = wh.findGenomeLocations(555933597,927287928,positions);
    cout << "got " << count << " positions:";
    for (unsigned int i=0; i<count; i++) cout << " " << positions[i];
    cout << "\n";


    return 0;
}

#endif
