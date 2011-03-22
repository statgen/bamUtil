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

#ifndef _MAPPERSECOLORSPACE_H_
#define _MAPPERSECOLORSPACE_H_
#include "MapperSE.h"

class MapperSEColorSpace: public MapperSE
{
public:
    MapperSEColorSpace();
    void MapSingleRead();

    int test(int testNum, const char *read, const char *qual, char direction, int chr, genomeIndex_t index, int misMatches, int quality);

    // the following two are not claimed to be pure virtual,
    // since MapperSEColorSpace class has not ready to implement these 2 functions yet.
    void MapSingleReadGapped();
    void MapSingleReadUnGapped();
    bool evalSEColorSpace(ReadIndexer &indexer,
                          genomeIndex_t genomeMatchPosition,
                          unsigned int whichWord);
};


#endif /* _MAPPERSECOLORSPACE_H_ */
