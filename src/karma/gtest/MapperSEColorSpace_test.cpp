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

class MapperSEColorSpaceTest: public ::testing::Test
{
public:
    GenomeSequence gs;
    GenomeSequence csgs;
    WordIndex wi;
    MapperSE* seMapper;

protected:
    virtual void SetUp()
    {
        if (gs.open(false) || csgs.open(true) || wi.open(csgs))
            error_at_line(1, 1, __FILE__,  __LINE__, "Failed to SetUp() !");

        seMapper = new MapperSEColorSpace;
        seMapper->setGenomeSequence(&csgs);
        seMapper->setWordIndex(&wi);

    };
    virtual void TearDown()
    {
        if (!seMapper) delete seMapper;
    }
};


// Test color space, single end, exact match
TEST_F(MapperSEColorSpaceTest, CS_SE_exactMatch)
{
    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment =         "C21202221110033303332010301313022211";
    String expectedReadFragment =  "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC";
    String dataQuality =          "!55555555555555555555555555555555566";
    String inputCigar="0M";

    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, readFragment, dataQuality) == 0);
    seMapper->MapSingleRead();
    std::ostringstream output;
    (seMapper->getBestMatch()).printColorSpace(output, &gs, &csgs, 0,
                                               readFragment, dataQuality,
                                               fragmentTag, true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, cszTag, csqTag;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual>> cszTag >> csqTag;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("0", flag.c_str());
    EXPECT_STREQ("7", S(chr));
    EXPECT_STREQ("110805955", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ(S(expectedReadFragment) , S(seq));
    // ignore this for now.
    // EXPECT_STREQ(S(MatchedReadBase::convertCSQuality(dataQuality)) , S(qual));
    String expectedCSZTag= "CS:Z:"+readFragment;
    String expectedCSQTag= "CQ:Z:"+dataQuality;
    EXPECT_STREQ(S(expectedCSZTag) , S(cszTag));
    EXPECT_STREQ(S(expectedCSQTag) , S(csqTag));

}

// // Test color space, single end, reverse of exact match
TEST_F(MapperSEColorSpaceTest, CS_SE_reverse)
{
    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment =                 "C21202221110033303332010301313022211";
    String expectedReadFragment =          "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC";
    //    String expectedReadFragment =          "TGAGAATGCAATTGGATATTATAAACACTCTTCA";
    String inverseComplementReadFragment = "T11122203131030102333033300111222021";
    String inverseReadFragment =           "A11122203131030102333033300111222021";
    String dataQuality =                   "!55555555555555555555555555555555566";
    String inverseDataQuality =             dataQuality.RightToLeft();
    String inputCigar="0M";

    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, inverseReadFragment, inverseDataQuality) == 0);
    seMapper->MapSingleRead();

    std::ostringstream output;
    (seMapper->getBestMatch()).printColorSpace(output, &gs, &csgs, 0,
                                               inverseReadFragment, inverseDataQuality,
                                               fragmentTag, true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, cszTag, csqTag;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual>> cszTag >> csqTag;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("16", flag.c_str()); // flag 0x0010 means reverse strand
    EXPECT_STREQ("7", S(chr));
    EXPECT_STREQ("110805955", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));

    // String expectedDataQuality = MatchedReadBase::convertCSQuality(dataQuality).RightToLeft();
    // EXPECT_STREQ(S(expectedDataQuality) , S(qual));
    String expectedCSZTag= "CS:Z:"+inverseReadFragment;
    String expectedCSQTag= "CQ:Z:"+inverseDataQuality;
    EXPECT_STREQ(S(expectedCSZTag) , S(cszTag));
    EXPECT_STREQ(S(expectedCSQTag) , S(csqTag));

    // test reverse complement
    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, inverseComplementReadFragment, inverseDataQuality) == 0);
    seMapper->MapSingleRead();

    output.seekp(0);
    (seMapper->getBestMatch()).printColorSpace(output, &gs, &csgs, 0,
                                               inverseComplementReadFragment, inverseDataQuality,
                                               fragmentTag, true, inputCigar);
    in.str(output.str());
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual>> cszTag >> csqTag;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("16", flag.c_str()); // flag 0x0010 means reverse strand
    EXPECT_STREQ("7", S(chr));
    EXPECT_STREQ("110805955", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ(S(expectedReadFragment) , S(seq));
    //EXPECT_STREQ(S(expectedDataQuality) , S(qual));
    expectedCSZTag= "CS:Z:"+inverseComplementReadFragment;
    expectedCSQTag= "CQ:Z:"+inverseDataQuality;
    EXPECT_STREQ(S(expectedCSZTag) , S(cszTag));
    EXPECT_STREQ(S(expectedCSQTag) , S(csqTag));

}

// Test color space, single end, 1 mismatch
TEST_F(MapperSEColorSpaceTest, CS_SE_1MisMatch)
{
    // now only for single end

    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment =        "C21202221110033303332010301313022211";
    String expectedReadFragment = "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC"; //len=35
    String dataQuality =         "!55555555555555555555555555555555555";
    String inputCigar="0M";

    String mutatedReadFragment;
    String mutatedReadFragmentWithPrimer;
    String mutatedCSReadFragment;
    char BASES[]={'A','T','G','C','N'};

    // always notice that color space with primer is 1 bp longer than color space read
    for (int i = 0; i< readFragment.Length() - 1 ; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            mutatedReadFragment = expectedReadFragment;
            if (mutatedReadFragment[i]==BASES[j])
                continue;
            mutatedReadFragment[i]=BASES[j];
            mutatedReadFragmentWithPrimer="A" + mutatedReadFragment;
            mutatedCSReadFragment=MatchedReadBase::convertBSRead(mutatedReadFragmentWithPrimer);
            ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, mutatedCSReadFragment, dataQuality) == 0);
            seMapper->MapSingleRead();

            std::ostringstream output;
            (seMapper->getBestMatch()).printColorSpace(output, &gs, &csgs, 0,
                                                       mutatedCSReadFragment, dataQuality,
                                                       fragmentTag, true, inputCigar);

            std::istringstream in(output.str(),std::istringstream::in);
            std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, cszTag, csqTag;
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
            EXPECT_STRCASEEQ(S(mutatedReadFragment) , S(seq)) << "seq= " <<seq ;
            EXPECT_STREQ(S(MatchedReadBase::convertCSQuality(dataQuality)) , S(qual));

            do
            {
                in>> cszTag;
            }
            while (strncmp(cszTag.c_str(),"CS:Z:",5) && in);
            in >> csqTag;
            String expectedCSZTag= "CS:Z:"+mutatedCSReadFragment;
            String expectedCSQTag= "CQ:Z:"+dataQuality;
            EXPECT_STREQ(S(expectedCSZTag) , S(cszTag));
            EXPECT_STREQ(S(expectedCSQTag) , S(csqTag));

        }
    }
}

// Test color space, single end, non mappable
TEST_F(MapperSEColorSpaceTest, CS_SE_nonMappable)
{
    // now only for single end
    String fragmentTag  = "@fake_example";
    String readFragment =         "A11111000002222333331111222333330000";
    String expectedReadFragment =  "CACACCCCCCTCTCGCGCGTGTGAGATATATTTTT";
    String dataQuality =          "!55555555555555555555555555555555555";
    String inputCigar="0M";

    ASSERT_TRUE(seMapper->processReadAndQuality(fragmentTag, readFragment, dataQuality) == 0);
    seMapper->MapSingleRead();

    std::ostringstream output;
    (seMapper->getBestMatch()).printColorSpace(output, &gs, &csgs, 0,
                                               readFragment, dataQuality,
                                               fragmentTag, true, inputCigar);
    // printf("%s\n", output.str().c_str());
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, haTag, nbTag, erTag;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual >> haTag >> nbTag >> erTag;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("4", flag.c_str());
    EXPECT_STREQ("1", S(chr));
    EXPECT_STREQ("0", S(pos));
    EXPECT_STREQ("0" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ(S(expectedReadFragment) , S(seq));
    EXPECT_STREQ(S(MatchedReadBase::convertCSQuality(dataQuality)) , S(qual));
}

// TODO
// Test color space, single end, 1 insertion
// Test color space, single end, 1 deletion

#endif
