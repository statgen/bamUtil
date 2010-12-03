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

TEST(GenomeSequenceClass, forBaseSpace)
{
    GenomeSequence gs;
    gs.setReferenceName("../testdata/phiX.fa");
    gs.open(false);
    // test phiX
    // >1 phiX: http://www.genome.jp/dbget-bin/www_bget?refseq+NC_001422
    // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAA
    // ....
    // AATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA (last line)
    EXPECT_EQ('G',gs[0]);
    EXPECT_EQ('A',gs[1]);
    EXPECT_EQ('G',gs[2]);
    EXPECT_EQ('T',gs[3]);
    EXPECT_EQ('T',gs[4]);

    int total = gs.getNumberBases();
    EXPECT_EQ('A',gs[--total]);
    EXPECT_EQ('C',gs[--total]);
    EXPECT_EQ('G',gs[--total]);

    // test total chromosome numbers
    EXPECT_EQ(1, gs.getChromosomeCount());

    // test getChromosome
    EXPECT_EQ(0, gs.getChromosome("1"));
    total = gs.getNumberBases();
    EXPECT_EQ(5386, total);
    EXPECT_EQ((genomeIndex_t) 5386, gs.getChromosomeSize(0)) <<"Chromosome 1 Length";
}

TEST(GenomeSequenceClass, forColorSpace)
{
    GenomeSequence gs;
    gs.setReferenceName("../testdata/phiX.fa");
    gs.open(true); // open color space
    // test phiX
    // >1 phiX: http://www.genome.jp/dbget-bin/www_bget?refseq+NC_001422
    // GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAA
    // => bs2cs
    // G22100033233202013121331220210301112002302333002212312212320
    // ....
    // AATTTATGCGCGCTTCGATAAAAATGATTGGCGTATCCAACCTGCA (last line)
    // => bs2cs
    // A030033133333202323300003123010331332010102131

    EXPECT_EQ('N',gs[0]);
    EXPECT_EQ('2',gs[1]);
    EXPECT_EQ('2',gs[2]);
    EXPECT_EQ('1',gs[3]);
    EXPECT_EQ('0',gs[4]);
    EXPECT_EQ('0',gs[5]);

    int total = gs.getNumberBases();
    EXPECT_EQ(5386, total);
    EXPECT_EQ('1', gs[--total]);
    EXPECT_EQ('3', gs[--total]);
    EXPECT_EQ('1', gs[--total]);
    EXPECT_EQ('2', gs[--total]);
}

class MapperBaseTest: public ::testing::Test
{
public:
    GenomeSequence gs;
    GenomeSequence csgs;
    MatchedReadBase matchedRead;

protected:
    virtual void SetUp()
    {
        gs.setReferenceName("../testdata/phiX.fa");
        gs.open(false);
        csgs.setReferenceName("../testdata/phiX.fa");
        csgs.open(true);
    }
    virtual void TearDown()
    {
    }
};

TEST_F(MapperBaseTest, calibrateSequenceTest)
{
    std::string read_fragment = "";
    std::string data_quality = "";
    // in base space, read equals   ATGTTAACTTCTGCGTCATGGAAGCGATAAAACTCT (origin)
    //                              TACAATTGAAGACGCAGTACCTTCGCTATTTTGAGA (reverse)
    //                              AGAGTTTTATCGCTTCCATGACGCAGAAGTTAACAT (reverse complement)
    std::string cs_read_fragment = "A31103012022133121310202332330001222";
    std::string cs_data_quality = "!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";
    genomeIndex_t genomeMatchPosition = 0;
    bool isForwardStrand = false;

    matchedRead.calibrateSequence(read_fragment, data_quality,
                                  cs_read_fragment, cs_data_quality,
                                  &gs, &csgs,
                                  genomeMatchPosition,
                                  isForwardStrand);
    // EXPECT_EQ("GAGTTTTATCGCTTCCATGACGCAGAAGTTAACAT", read_fragment);
    //        AGAGTTTTATCGCTTCCATGACGCAGAAGTTAACAT (reverse complement of the read)
    //       GAGTTTTATCGCTTCCATGACGCAGAAGTTAACACTTTCGGATATTTCTGATGAGTCGAA (first line of phiX.fa)
    // unfinished
};

TEST_F(MapperBaseTest, fixBaseRangeTest)
{
    int start = 1;
    int end = 2;
    std::string read_fragment = "";
    std::string data_quality = "";
    std::string cs_read_fragment = "";
    std::string cs_data_quality = "";
    genomeIndex_t genomeMatchPosition = 0;
    bool isForwardStrand = true;

    matchedRead.fixBaseRange(start,  end, /// inclusive boundaries
                             read_fragment, data_quality,
                             cs_read_fragment, cs_data_quality,
                             &gs, &csgs,
                             genomeMatchPosition,
                             isForwardStrand);
    // EXPECT_EQ("" , read_fragment);
    // EXPECT_EQ("" , data_quality);
    // TODO: unfisihed
    SUCCEED();
};

TEST_F(MapperBaseTest, ColorSpaceRelatedFunction)
{
    // convertCSReadTest
    std::string input;
    std::string ret;
    input = "A00";
    EXPECT_EQ("AA", matchedRead.convertCSRead(input));
    EXPECT_EQ("AAA", matchedRead.convertCSRead(input, true));

    input = "A011223";
    EXPECT_EQ("AACAGAT", matchedRead.convertCSRead(input, true));

    input =   "G310022100002020303200022112.1.21.2.";
    EXPECT_EQ("GCAAAGACCCCCTTCCGGCTTTTCTGTCNNNNNNNN", matchedRead.convertCSRead(input,true));

    // test reverse, complement, reverse complement
    input = "ACCGT";
    EXPECT_EQ("TGCCA", matchedRead.getReverse(input));
    EXPECT_EQ("TGGCA", matchedRead.getComplement(input));
    EXPECT_EQ("ACGGT", matchedRead.getReverseComplement(input));

    input = "A011223";
    EXPECT_EQ("322110A", matchedRead.getReverse(input));
    EXPECT_EQ("T011223", matchedRead.getComplement(input));
    EXPECT_EQ("322110T", matchedRead.getReverseComplement(input));


};

//////////////////////////////////////////////////////////////////////
// Old code
//////////////////////////////////////////////////////////////////////
#if 0
class BaseSpaceTest: public ::testing::Test
{
public:
    GenomeSequence gs;
    WordIndex wi;
    WordHash whLeft;
    WordHash whRight;
    MapperSE* seMapper;
    MapperPE* peMapper;

protected:
    virtual void SetUp()
    {
        if (gs.open(false) || wi.open(gs))
            error(1, 1,  "Failed to SetUp at __LINE__!");
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
        if (!peMapper) delete peMapper;
    }
};

TEST_F(BaseSpaceTest, testGetReadAndQuality)
{
    String filename="/home/zhanxw/compareMapSoft/testSuite/singleBS.fastq";
    IFILE f = ifopen(filename, "rb");
    EXPECT_EQ(true, f !=NULL) << "Reads file [ "<< filename.c_str()<< "can not be opened\n";

    seMapper->getReadAndQuality(f);
    EXPECT_STREQ("TTTTTATTTATATATATATAAAATATTGACAAAGA", (seMapper->forward).read.c_str());
    EXPECT_STREQ("TCTTTGTCAATATTTTATATATATATAAATAAAAA", (seMapper->backward).read.c_str());
    ifclose(f);
}


TEST_F(BaseSpaceTest, printBaseSpace)
{
    // now only for single end
    std::ostringstream output;
    String filename="/home/zhanxw/compareMapSoft/testSuite/singleBS.fastq";
    IFILE f = ifopen(filename, "rb");
    ASSERT_TRUE(seMapper->getReadAndQuality(f) == 0);
    seMapper->MapSingleRead();

    String fragmentTag="@Chromosome_4__132622643_Genome_0822325337";
    String inputCigar="0M";
    (seMapper->getBestMatch()).print(output, 0, fragmentTag, true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
    EXPECT_STREQ("Chromosome_4__132622643_Genome_0822325337", S(tag));
    EXPECT_STREQ("0", S(flag));
    EXPECT_STREQ("4", S(chr));
    EXPECT_STREQ("132622643", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ("TTTTTATTTATATATATATAAAATATTGACAAAGA" , S(seq));
    EXPECT_STREQ("55555555555555555555555555555555555" , S(qual));

    ifclose(f);
}

// Test base space, paired end, exact match
TEST_F(BaseSpaceTest, BS_SE_exactMatch)
{
    String fragmentTag  = "@Chromosome_7__110805955_Genome_1343539570";
    String readFragment = "TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC";
    String dataQuality =  "55555555555555555555555555555555555";
    String inputCigar="0M";

    ASSERT_TRUE(peMapper->processReadAndQuality(fragmentTag, readFragment, dataQuality) == 0);

    seMapper->MapSingleRead();

    std::ostringstream output;
    (seMapper->getBestMatch()).print(output, 0, fragmentTag, true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("0", flag.c_str());
    EXPECT_TRUE("7"== chr);
    EXPECT_TRUE("110805955"== pos);
    EXPECT_TRUE("100" ==  mapq);
    EXPECT_TRUE("0M" ==  cigar);
    EXPECT_TRUE("=" ==  mchr);
    EXPECT_TRUE("0" ==  mpos);
    EXPECT_TRUE("0" ==  isize);
    EXPECT_TRUE("TGAAGAGTGTTTATAATATCCAATTGCATTCTCAC" == seq);
    EXPECT_TRUE("55555555555555555555555555555555555" == qual);

}

// Test base space, paired end, 1 mismatch
// Test base space, paired end, non mappable
// Test base space, paired end, 1 insertion
// Test base space, paired end, 1 deletion



class ColorSpaceTest: public ::testing::Test
{
public:
    GenomeSequence gs;
    GenomeSequence csgs;
    WordIndex cswi;

protected:
    virtual void SetUp()
    {
        if (gs.open(false) || csgs.open(true) || cswi.open(csgs))
            error(1, 1,  "Failed to SetUp at __LINE__!");
    };
};

TEST_F(ColorSpaceTest, testGetReadAndQuality)
{
    String filename="/home/zhanxw/compareMapSoft/testSuite/singleCS.fastq";
    IFILE f = ifopen(filename, "rb");

    EXPECT_EQ(true, f !=NULL) << "Reads file [ "<< filename.c_str()<< "can not be opened\n";

    MapperSE* mapper=new MapperSEColorSpace;
    GenomeSequence gs;
    gs.open(true);
    WordIndex wi;
    wi.open(gs);
    mapper->setGenomeSequence(&gs);
    mapper->setWordIndex(&wi);
    mapper->getReadAndQuality(f);
    EXPECT_STREQ("AAAATTAATTTTTTTTTTTAAATTTACGCCAAGG", (mapper->forward).read.c_str());
    EXPECT_STREQ("GGAACCGCATTTAAATTTTTTTTTTTAATTAAAA", (mapper->backward).read.c_str());

    ifclose(f);
}


TEST_F(ColorSpaceTest, printColorSpace)
{
    // now only for single end
    std::ostringstream output;
    MapperSE* mapper=new MapperSEColorSpace;
    mapper->setGenomeSequence(&csgs);
    mapper->setWordIndex(&cswi);

    String filename="/home/zhanxw/compareMapSoft/testSuite/singleCS.fastq";
    IFILE f = ifopen(filename, "rb");
    mapper->getReadAndQuality(f);
    mapper->MapSingleRead();

    String fragmentTag="@Chromosome_4__132622643_Genome_0822325337";
    String inputCigar="0M";
    (mapper->getBestMatch()).printColorSpace(output, &gs, &csgs, 0,
            mapper->originalCSRead,
            mapper->originalCSQual,
            fragmentTag,
            true, inputCigar);
    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, hatag, cstag, cqtag;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual>>hatag>>cstag>>cqtag;
    EXPECT_STREQ("Chromosome_4__132622643_Genome_0822325337", S(tag));
    EXPECT_STREQ("0", S(flag));
    EXPECT_STREQ("4", S(chr));
    EXPECT_STREQ("132622643", S(pos));
    EXPECT_STREQ("90" , S(mapq)) << "mapq = "<<mapq;
    EXPECT_STREQ("0M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("0" , S(mpos));
    EXPECT_STREQ("0" , S(isize));
    EXPECT_STREQ("TTTTTATTTATATATATATAAAATATTGACAAAGA" , S(seq));
    EXPECT_STREQ("!1@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!" , S(qual));
    EXPECT_STREQ("HA:i:890" , S(hatag));
    EXPECT_STREQ("CS:Z:A30000330033333333333000333012110022", S(cstag));
    EXPECT_STREQ("CQ:Z:!1111111111111111111111111111111111", S(cqtag));

    delete mapper;
    ifclose(f);
}

// Test colorspace, single end, exact match
// Test colorspace, single end, 1 mismatch
// Test colorspace, single end, non mappable
// Test colorspace, single end, 1 insertion
// Test colorspace, single end, 1 deletion

// Test colorspace, paired end, exact match
// Test colorspace, paired end, 1 mismatch
// Test colorspace, paired end, non mappable
// Test colorspace, paired end, 1 insertion
// Test colorspace, paired end, 1 deletion
#endif
