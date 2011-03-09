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

#include "GenomeSequence.h"
#include "MapperBase.h"
#include "MapperSEBaseSpace.h"
#include "MapperSEColorSpace.h"

#include <error.h>
#include <gtest/gtest.h>


TEST(IntegrityCheck, forBaseSpace)
{
    GenomeSequence gs;
    gs.setReferenceName("../testdata/phiX.fa");
    if (gs.open(false))
        error(1,1, "Cannot open reference genome at __LINE__");
    genomeIndex_t index=0;

    FILE* f=fopen("../testdata/phiX.fa","rt");
    char data;
    while (!feof(f))
    {
        fscanf(f,"%c",&data);
        switch (data)
        {
        case '\n':
        case '\r':
            continue;
        break;
        case '>':
            while (data!='\n' && data!='\r')
                fscanf(f,"%c",&data);
            break;
        default:
            ASSERT_EQ(GenomeSequence::base2int[(int)gs[index]],
                      GenomeSequence::base2int[(int)data])
                << "gs[index]=" << gs[index]
                << "data=" << data;
            index++;
        }
    }
    fclose(f);
}

TEST(IntegrityCheck, forColorSpace)
{
    // base codes: 0-> A,    1-> C,     2-> G,      3-> T
    // colorspace: 0-> blue, 1-> green, 2-> oragne, 3->red
    //
    const char fromBase2CS[] =
        {
            /* 0000 */ 0,   // A->A
            /* 0001 */ 1,   // A->C
            /* 0010 */ 2,   // A->G
            /* 0011 */ 3,   // A->T
            /* 0100 */ 1,   // C->A
            /* 0101 */ 0,   // C->C
            /* 0110 */ 3,   // C->G
            /* 0111 */ 2,   // C->T
            /* 1000 */ 2,   // G->A
            /* 1001 */ 3,   // G->C
            /* 1010 */ 0,   // G->G
            /* 1011 */ 1,   // G->T
            /* 1100 */ 3,   // T->A
            /* 1101 */ 2,   // T->C
            /* 1110 */ 1,   // T->G
            /* 1111 */ 0,   // T->T
        };

    GenomeSequence gs;
    gs.setReferenceName("../testdata/phiX.fa");
    if (gs.open(true))
        error(1,1, "Cannot open reference genome at __LINE__");
    genomeIndex_t index=0;

    FILE* f=fopen("../testdata/phiX.fa","rt");
    char data;
    signed char lastBase = 5;    // for converting the reference to colorspace, the first bp is always 5 (in base space it is 'N')
    while (!feof(f))
    {
        fscanf(f,"%c",&data);
        switch (data)
        {
        case '\n':
        case '\r':
            continue;
        break;
        case '>':
            while (data!='\n' && data!='\r')
                fscanf(f,"%c",&data);
            break;
        default:
            // so we don't write a colorspace value when we
            // get the first base value.
            //
            // On second and subsequent bases, write based on
            // the index table
            //
            char thisBase = GenomeSequence::base2int[(int)(data)];
            if (lastBase>=0)
            {
                char color;
                if (lastBase>3 || thisBase>3) color=4;
                else color = fromBase2CS[(int)(lastBase<<2 | thisBase)];
                // re-use the int to base, because ::set expects a base char (ATCG), not
                // a color code (1234).  It should only matter on final output.
                ASSERT_EQ(GenomeSequence::base2int[(int)gs[index]],
                          color)
                    << "gs[index]=" << gs[index]
                    << "data=" << data
                    << "index=" <<index;
            }
            lastBase = thisBase;

            index++;
        }
    }
    fclose(f);
}

TEST(IntegrityCheck_WordIndex, forBaseSpace)
{
#if 0
    GenomeSequence gs;
    WordIndex wi;
    if (gs.open(false) || wi.open(gs))
        error(1,1, "Cannot open reference genome at __LINE__");

    std::cout<<"If you see anything below, then word index is probably corrupt!"<<std::endl;
    wi.testHashIndices(gs);
    std::cout<<"--------------- END ---------------"<<std::endl;
#endif
    SUCCEED();
}

TEST(IntegrityCheck_WordIndex, forColorSpace)
{
#if 0
    GenomeSequence gs;
    WordIndex wi;
    if (gs.open(true) || wi.open(gs))
        error(1,1, "Cannot open reference genome at __LINE__");

    std::cout<<"If you see anything below, then word index is probably corrupt!"<<std::endl;
    wi.testHashIndices(gs);
    std::cout<<"--------------- END ---------------"<<std::endl;
#endif
    SUCCEED();
}

