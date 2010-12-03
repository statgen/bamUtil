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
#include "MapperSE.h"
#include "MapperSEBaseSpace.h"
#include "MapperSEColorSpace.h"
#include "MapperPE.h"
#include "MapperPEBaseSpace.h"
#include "MapperPEColorSpace.h"

#include <error.h>
#include <gtest/gtest.h>

#if 0
#define S(x) ((x).c_str())

class MapperSEBaseSpaceTest: public ::testing::Test
{
public:
    GenomeSequence gs;
    WordIndex wi;
    WordHash whLeft;
    WordHash whRight;
    MapperSE* seMapper;

protected:
    virtual void SetUp()
    {
        if (gs.open(false) || wi.open(gs))
            error_at_line(1, 1, __FILE__,  __LINE__, "Failed to SetUp() !");
        if (whLeft.open(gs.getBaseFilename() + ".umwhl"))
        {
            std::cerr << "failed to open left word hash " << gs.getBaseFilename() + ".umwhl"  << ".\n";
            exit(1);
        }
        if (whRight.open(gs.getBaseFilename() + ".umwhr"))
        {
            std::cerr << "failed to open right word hash " << gs.getBaseFilename() + ".umwhl"  << ".\n";
            exit(1);
        }

        seMapper = new MapperSEBaseSpace;

        seMapper->setGenomeSequence(&gs);
        seMapper->setWordIndex(&wi);
        seMapper->setWordHashLeft(&whLeft);
        seMapper->setWordHashRight(&whRight);

    };
    virtual void TearDown()
    {
        if (!seMapper) delete seMapper;
    }
};


// Test base space, single end, exact match
TEST_F(MapperSEBaseSpaceTest, BS_SE_exactMatch)
{
    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment = "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC";
    String dataQuality =  "55555555555555555555555555555555556";
    String inputCigar="0M";

    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, readFragment, dataQuality) == 0);
    seMapper->MapSingleRead();

    std::ostringstream output;
    (seMapper->getBestMatch()).print(output, 0, fragmentTag, true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("0", flag.c_str());
    EXPECT_STREQ("7", S(chr));
    EXPECT_STREQ("110805955", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ(S(readFragment) , S(seq));
    EXPECT_STREQ(S(dataQuality) , S(qual));

}
// Test base space, single end, reverse of exact match
TEST_F(MapperSEBaseSpaceTest, BS_SE_reverse)
{
    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment = "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC";
    String dataQuality =  "55555555555555555555555555555555556";
    String inverseReadFragment = readFragment.RightToLeft();
    String inverseDataQuality = dataQuality.RightToLeft();
    String inputCigar="0M";

    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, inverseReadFragment, inverseDataQuality) == 0);
    seMapper->MapSingleRead();

    std::ostringstream output;
    (seMapper->getBestMatch()).print(output, 0, fragmentTag, true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("4", flag.c_str()); // unmappable
    EXPECT_STREQ("1", S(chr));
    EXPECT_STREQ("1", S(pos));
    EXPECT_STREQ("0" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ(S(inverseReadFragment) , S(seq));
    EXPECT_STREQ(S(inverseDataQuality) , S(qual));

}

// Test base space, single end, 1 mismatch
TEST_F(MapperSEBaseSpaceTest, BS_SE_1MisMatch)
{
    // now only for single end

    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment = "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC";
    String dataQuality =  "55555555555555555555555555555555555";
    String inputCigar="0M";

    String mutatedReadFragment;
    char BASES[]={'A','T','G','C','N'};
    for (int i = 0; i< readFragment.Length() ; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mutatedReadFragment = readFragment;
            if (mutatedReadFragment[i]==BASES[j])
                continue;
            mutatedReadFragment[i]=BASES[j];
            ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, mutatedReadFragment, dataQuality) == 0);
            seMapper->MapSingleRead();

            std::ostringstream output;
            (seMapper->getBestMatch()).print(output, 0, fragmentTag, true, inputCigar);
            std::istringstream in(output.str(),std::istringstream::in);
            std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
            in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
            EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
            EXPECT_STREQ("0", flag.c_str());
            EXPECT_STREQ("7", S(chr));
            EXPECT_STREQ("110805955", S(pos));
            EXPECT_STREQ("100" , S(mapq));
            EXPECT_STREQ("0M" , S(cigar));
            EXPECT_STREQ("=" , S(mchr));
            EXPECT_STREQ("0" , S(mpos));
            EXPECT_STREQ("0" , S(isize));
            EXPECT_STRCASEEQ(mutatedReadFragment , S(seq)) << "seq= " <<seq;
            EXPECT_STREQ(S(dataQuality) , S(qual));
        }
    }
}

// Test base space, single end, non mappable
TEST_F(MapperSEBaseSpaceTest, BS_SE_nonMappable)
{
    // now only for single end
    String fragmentTag  = "@@Chr_6__053181553_i_0400/2";
    String readFragment = "CTGCACTCCAGCCTGGGTGACAGAGCGAGACTCTG";
    String dataQuality =  "55555555555555555555555555555555555";
    String inputCigar="0M";

    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, readFragment, dataQuality) == 0);
    seMapper->MapSingleRead();

    std::ostringstream output;
    (seMapper->getBestMatch()).print(output, 0, fragmentTag, true, inputCigar);
    // printf("%s\n", output.str().c_str());
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, haTag, nbTag, erTag;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual >> haTag >> nbTag >> erTag;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("4", flag.c_str());
    EXPECT_STREQ("1", S(chr));
    EXPECT_STREQ("1", S(pos));
    EXPECT_STREQ("0" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ(S(readFragment) , S(seq));
    EXPECT_STREQ(S(dataQuality) , S(qual));
    EXPECT_STREQ("HA:i:22", S(haTag));
    EXPECT_STREQ("NB:i:6", S(nbTag));
    EXPECT_STREQ("ER:Z:duplicates", S(erTag));


}

// TODO
// Test base space, single end, 1 insertion
// Test base space, single end, 1 deletion

#endif

