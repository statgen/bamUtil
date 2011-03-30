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

class MapperPEColorSpaceTest: public ::testing::Test
{
public:
    GenomeSequence gs;
    GenomeSequence csgs;
    WordIndex wi;
    WordHash whLeft;
    WordHash whRight;
    MapperPE* peMapperA;
    MapperPE* peMapperB;

    String fragmentTag ;
    String fragmentTag2 ;

    String readFragment ;
    String readFragmentBS ;
    String expectedReadFragmentBS;
    String readFragment2 ;
    String readFragment2BS ;
    String expectedReadFragment2BS;

    String inverseComplementReadFragment ;
    String inverseComplementReadFragment2 ;

    String dataQuality ;
    String inverseDataQuality ;
    String dataQuality2 ;
    String inverseDataQuality2 ;

protected:
    virtual void SetUp()
    {
        if (gs.open(false) || csgs.open(true) || wi.open(csgs))
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

        MapperUserOption mapperOptions;
        mapperOptions.showReferenceBases=true;

        peMapperA = new MapperPEColorSpace;
        peMapperA->setGenomeSequence(&csgs);
        peMapperA->setWordIndex(&wi);
        peMapperA->setWordHashLeft(&whLeft);
        peMapperA->setWordHashRight(&whRight);
        peMapperA->setMapperOptions(mapperOptions);

        peMapperB = new MapperPEColorSpace;
        peMapperB->setGenomeSequence(&csgs);
        peMapperB->setWordIndex(&wi);
        peMapperB->setWordHashLeft(&whLeft);
        peMapperB->setWordHashRight(&whRight);
        peMapperB->setMapperOptions(mapperOptions);


        fragmentTag  = "@Chr_X__113361308_i_0400/1";
        fragmentTag2  = "@Chr_X__113361308_i_0400/2";

        readFragment =            "A33113120333112111331211033201302132";
        readFragmentBS =           "TACATGAATATGTCACATACTGTTATCCATTCATC";
        expectedReadFragmentBS =  readFragmentBS;
        readFragment2 =           "A11023132212002203312212111033032121";
        readFragment2BS=           "CAAGCATCTGAAAGAATACTCAGTGTTATTAGTCA";
        expectedReadFragment2BS =  "TGACTAATAACACTGAGTATTCTTTCAGATGCTTG"; // the forward strand seen in ref genome

        dataQuality  =            "!35555555555555555555555555555555557";
        inverseDataQuality = dataQuality.RightToLeft();
        dataQuality2 =            "!45555555555555555555555555555555556";
        inverseDataQuality2 = dataQuality2.RightToLeft();

    };
    virtual void TearDown()
    {
        if (!peMapperA) delete peMapperA;
        if (!peMapperB) delete peMapperB;
    }
};


// Test color space, paired end, exact match
TEST_F(MapperPEColorSpaceTest, CS_PE_exactMatch)
{

    ASSERT_TRUE(peMapperA->processReadAndQuality(fragmentTag,  readFragment,  dataQuality) == 0);
    ASSERT_TRUE(peMapperB->processReadAndQuality(fragmentTag2, readFragment2, dataQuality2) == 0);
    peMapperA->mapReads(peMapperB);

    std::ostringstream output;
    peMapperA->printCSBestReads(output, &gs, &csgs, peMapperB);

    std::istringstream in(output.str(),std::istringstream::in);
    std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual, cszTag, csqTag;
    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual>> cszTag >> csqTag;
    EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
    EXPECT_STREQ("33", flag.c_str());
    EXPECT_STREQ("X", S(chr));
    EXPECT_STREQ("113361308", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("34M" , S(cigar)); // TODO: this will be improved
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("113361708" , S(mpos));
    EXPECT_STREQ("-435" , S(isize));
    EXPECT_STREQ(S(expectedReadFragmentBS) , S(seq));
    EXPECT_STREQ(S(MatchedReadBase::convertCSQuality(dataQuality)) , S(qual));
    String expectedCSZTag= "CS:Z:"+readFragment;
    String expectedCSQTag= "CQ:Z:"+dataQuality;
    EXPECT_STREQ(S(expectedCSZTag) , S(cszTag));
    EXPECT_STREQ(S(expectedCSQTag) , S(csqTag));

    in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual>> cszTag >> csqTag;
    EXPECT_STREQ(fragmentTag2.c_str()+1, tag.c_str());
    EXPECT_STREQ("17", flag.c_str());
    EXPECT_STREQ("X", S(chr));
    EXPECT_STREQ("113361708", S(pos));
    EXPECT_STREQ("100" , S(mapq));
    EXPECT_STREQ("34M" , S(cigar));
    EXPECT_STREQ("=" , S(mchr));
    EXPECT_STREQ("113361308" , S(mpos));
    EXPECT_STREQ("435" , S(isize));
    EXPECT_STREQ(S(expectedReadFragment2BS), S(seq));
    EXPECT_STREQ(S(MatchedReadBase::convertCSQuality(inverseDataQuality2)) , S(qual));

    expectedCSZTag= "CS:Z:"+readFragment2;
    expectedCSQTag= "CQ:Z:"+dataQuality2;
    EXPECT_STREQ(S(expectedCSZTag) , S(cszTag));
    EXPECT_STREQ(S(expectedCSQTag) , S(csqTag));

}

// // Test color space, paired end, reverse of exact match
// // Suppose input reads are readA, readB, we will test mapping :
// //  1. reverse complement of readA, reverse complement of readB
// //  2. readB, readA
// //
// TEST_F(MapperPEColorSpaceTest, CS_PE_reverse)
// {
//     // Test case 1: reverse complement of readA, reverse complement of readB
//     String fragmentTag  = "@Chr_X__113361308_i_0400/1";
//     String readFragment = "TACATGAATATGTCACATACTGTTATCCATTCATC";
//     String inverseComplementReadFragment = "GATGAATGGATAACAGTATGTGACATATTCATGTA";
//     String fragmentTag2  = "@Chr_X__113361308_i_0400/2";
//     String readFragment2 = "CAAGCATCTGAAAGAATACTCAGTGTTATTAGTCA";
//     String inverseComplementReadFragment2 = "TGACTAATAACACTGAGTATTCTTTCAGATGCTTG";
//     String dataQuality  =  "35555555555555555555555555555555557";
//     String inverseDataQuality = dataQuality.RightToLeft();
//     String dataQuality2 =  "45555555555555555555555555555555556";
//     String inverseDataQuality2 = dataQuality2.RightToLeft();

//     ASSERT_TRUE( peMapperA->processReadAndQuality(fragmentTag,  inverseComplementReadFragment,  dataQuality) == 0);
//     ASSERT_TRUE( peMapperB->processReadAndQuality(fragmentTag2, inverseComplementReadFragment2, dataQuality2) == 0);
//     peMapperA->mapReads(peMapperB);

//     std::ostringstream output;
//     peMapperA->printCSBestReads(output, peMapperB);

//     std::istringstream in(output.str(),std::istringstream::in);
//     std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
//     in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//     EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
//     EXPECT_STREQ("17", flag.c_str());
//     EXPECT_STREQ("X", S(chr));
//     EXPECT_STREQ("113361308", S(pos));
//     EXPECT_STREQ("100" , S(mapq));
//     EXPECT_STREQ("35M" , S(cigar));
//     EXPECT_STREQ("=" , S(mchr));
//     EXPECT_STREQ("113361708" , S(mpos));
//     EXPECT_STREQ("-435" , S(isize));
//     EXPECT_STREQ(S(readFragment) , S(seq));
//     EXPECT_STREQ(S(inverseDataQuality) , S(qual));

//     in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//     EXPECT_STREQ(fragmentTag2.c_str()+1, tag.c_str());
//     EXPECT_STREQ("33", flag.c_str());
//     EXPECT_STREQ("X", S(chr));
//     EXPECT_STREQ("113361708", S(pos));
//     EXPECT_STREQ("100" , S(mapq));
//     EXPECT_STREQ("35M" , S(cigar));
//     EXPECT_STREQ("=" , S(mchr));
//     EXPECT_STREQ("113361308" , S(mpos));
//     EXPECT_STREQ("435" , S(isize));
//     EXPECT_STREQ(S(inverseComplementReadFragment2), S(seq));
//     EXPECT_STREQ(S(dataQuality2) , S(qual));

//     // Test case 2: readB, readA
//     ASSERT_TRUE( peMapperA->processReadAndQuality(fragmentTag,  readFragment,  dataQuality) == 0);
//     ASSERT_TRUE( peMapperB->processReadAndQuality(fragmentTag2, readFragment2, dataQuality2) == 0);
//     peMapperB->mapReads(peMapperA);

//     output.seekp(0);
//     peMapperB->printCSBestReads(output, peMapperA);

//     in.str(output.str());
//     in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//     EXPECT_STREQ(fragmentTag2.c_str()+1, tag.c_str());
//     EXPECT_STREQ("17", flag.c_str());
//     EXPECT_STREQ("X", S(chr));
//     EXPECT_STREQ("113361708", S(pos));
//     EXPECT_STREQ("100" , S(mapq));
//     EXPECT_STREQ("35M" , S(cigar));
//     EXPECT_STREQ("=" , S(mchr));
//     EXPECT_STREQ("113361308" , S(mpos));
//     EXPECT_STREQ("435" , S(isize));
//     EXPECT_STREQ(S(inverseComplementReadFragment2), S(seq));
//     EXPECT_STREQ(S(inverseDataQuality2) , S(qual));

//     in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//     EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
//     EXPECT_STREQ("33", flag.c_str());
//     EXPECT_STREQ("X", S(chr));
//     EXPECT_STREQ("113361308", S(pos));
//     EXPECT_STREQ("100" , S(mapq));
//     EXPECT_STREQ("35M" , S(cigar));
//     EXPECT_STREQ("=" , S(mchr));
//     EXPECT_STREQ("113361708" , S(mpos));
//     EXPECT_STREQ("-435" , S(isize));
//     EXPECT_STREQ(S(readFragment) , S(seq));
//     EXPECT_STREQ(S(dataQuality) , S(qual));
// }

// // Test color space, paired end, 1 mismatch
// TEST_F(MapperPEColorSpaceTest, CS_PE_1MisMatch)
// {
//     // Mutate first read
//     String mutatedReadFragment;
//     char BASES[]={'A','T','G','C','N'};
//     for (int i = 0; i< readFragment.Length() ; i++)
//     {
//         for (int j = 0; j < 4; j++)
//         {
//             mutatedReadFragment = readFragment;
//             if (mutatedReadFragment[i]==BASES[j])
//                 continue;
//             mutatedReadFragment[i]=BASES[j];
//             ASSERT_TRUE( peMapperA->processReadAndQuality(fragmentTag,  mutatedReadFragment,  dataQuality) == 0);
//             ASSERT_TRUE( peMapperB->processReadAndQuality(fragmentTag2, readFragment2, dataQuality2) == 0);
//             peMapperA->mapReads(peMapperB);

//             std::ostringstream output;
//             peMapperA->printCSBestReads(output, peMapperB);

//             std::istringstream allInput(output.str(),std::istringstream::in);
//             std::istringstream in;
//             std::string line;
//             getline(allInput, line);
//             in.str(line);
//             std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
//             in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//             EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
//             EXPECT_STREQ("33", flag.c_str());
//             EXPECT_STREQ("X", S(chr));
//             EXPECT_STREQ("113361308", S(pos));
//             EXPECT_STREQ("100" , S(mapq));
//             EXPECT_STREQ("35M" , S(cigar));
//             EXPECT_STREQ("=" , S(mchr));
//             EXPECT_STREQ("113361708" , S(mpos));
//             EXPECT_STREQ("-435" , S(isize));
//             EXPECT_STRCASEEQ(S(mutatedReadFragment), S(seq));
//             EXPECT_STREQ(S(dataQuality) , S(qual));

//             getline(allInput, line);
//             in.str(line);
//             in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//             EXPECT_STREQ(fragmentTag2.c_str()+1, tag.c_str());
//             EXPECT_STREQ("17", flag.c_str());
//             EXPECT_STREQ("X", S(chr));
//             EXPECT_STREQ("113361708", S(pos));
//             EXPECT_STREQ("100" , S(mapq));
//             EXPECT_STREQ("35M" , S(cigar));
//             EXPECT_STREQ("=" , S(mchr));
//             EXPECT_STREQ("113361308" , S(mpos));
//             EXPECT_STREQ("435" , S(isize));
//             EXPECT_STREQ(S(inverseComplementReadFragment2), S(seq));
//             EXPECT_STREQ(S(inverseDataQuality2) , S(qual));
//         }
//     }

//     String invertComplementMutatedReadFrament;
//     for (int i = 0; i< readFragment2.Length() ; i++)
//     {
//         for (int j = 0; j < 4; j++)
//         {
//             mutatedReadFragment = readFragment2;
//             if (mutatedReadFragment[i]==BASES[j])
//                 continue;
//             mutatedReadFragment[i]=BASES[j];
//             invertComplementMutatedReadFrament=mutatedReadFragment.RightToLeft();
//             for (int temp=0; temp<invertComplementMutatedReadFrament.Length(); temp++)
//                 invertComplementMutatedReadFrament[temp]=gs.BasePair(invertComplementMutatedReadFrament[temp]);
//             ASSERT_TRUE( peMapperA->processReadAndQuality(fragmentTag,  readFragment,  dataQuality) == 0);
//             ASSERT_TRUE( peMapperB->processReadAndQuality(fragmentTag2, mutatedReadFragment, dataQuality2) == 0);
//             peMapperA->mapReads(peMapperB);

//             std::ostringstream output;
//             peMapperA->printCSBestReads(output, peMapperB);

//             std::istringstream allInput(output.str(),std::istringstream::in);
//             std::istringstream in;
//             std::string line;
//             getline(allInput, line);
//             in.str(line);
//             std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
//             in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//             EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
//             EXPECT_STREQ("33", flag.c_str());
//             EXPECT_STREQ("X", S(chr));
//             EXPECT_STREQ("113361308", S(pos));
//             EXPECT_STREQ("100" , S(mapq));
//             EXPECT_STREQ("35M" , S(cigar));
//             EXPECT_STREQ("=" , S(mchr));
//             EXPECT_STREQ("113361708" , S(mpos));
//             EXPECT_STREQ("-435" , S(isize));
//             EXPECT_STREQ(S(readFragment), S(seq));
//             EXPECT_STREQ(S(dataQuality) , S(qual));

//             getline(allInput, line);
//             in.str(line);
//             in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//             EXPECT_STREQ(fragmentTag2.c_str()+1, tag.c_str());
//             EXPECT_STREQ("17", flag.c_str());
//             EXPECT_STREQ("X", S(chr));
//             EXPECT_STREQ("113361708", S(pos));
//             EXPECT_STREQ("100" , S(mapq));
//             EXPECT_STREQ("35M" , S(cigar));
//             EXPECT_STREQ("=" , S(mchr));
//             EXPECT_STREQ("113361308" , S(mpos));
//             EXPECT_STREQ("435" , S(isize));
//             EXPECT_STRCASEEQ(S(invertComplementMutatedReadFrament), S(seq));
//             EXPECT_STREQ(S(inverseDataQuality2) , S(qual));
//         }
//     }
// }

// // Test color space, single end, non mappable
// TEST_F(MapperPEColorSpaceTest, CS_PE_nonMappable)
// {
//     String fragmentTag  = "@Chr_6__053181553_i_0400/1";
//     String readFragment = "CGCTGGATGAGATGAGAAAGGGCATTTATAAGTTC";
//     String fragmentTag2  = "@Chr_6__053181553_i_0400/2";
//     String readFragment2 = "CTGCACTCCAGCCTGGGTGACAGAGCGAGACTCTG";
//     String dataQuality =  "55555555555555555555555555555555555";

//     ASSERT_TRUE( peMapperA->processReadAndQuality(fragmentTag,  readFragment,  dataQuality) == 0);
//     ASSERT_TRUE( peMapperB->processReadAndQuality(fragmentTag2, readFragment2, dataQuality) == 0);
//     peMapperA->mapReads(peMapperB);

//     std::ostringstream output;
//     peMapperA->printCSBestReads(output, peMapperB);

//     std::istringstream allInput(output.str(),std::istringstream::in);
//     std::istringstream in;
//     std::string line;
//     getline(allInput, line);
//     in.str(line);
//     std::string tag, flag, chr, pos, mapq, cigar, mchr, mpos, isize, seq, qual;
//     in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//     EXPECT_STREQ(fragmentTag.c_str()+1, tag.c_str());
//     EXPECT_STREQ("13", flag.c_str());
//     EXPECT_STREQ("1", S(chr));
//     EXPECT_STREQ("1", S(pos));
//     EXPECT_STREQ("0" , S(mapq));
//     EXPECT_STREQ("35M" , S(cigar));
//     EXPECT_STREQ("=" , S(mchr));
//     EXPECT_STREQ("1" , S(mpos));
//     EXPECT_STREQ("35" , S(isize));
//     EXPECT_STREQ(S(readFragment) , S(seq));
//     EXPECT_STREQ(S(dataQuality) , S(qual));

//     getline(allInput, line);
//     in.str(line);
//     in>>tag>> flag>> chr>> pos>> mapq>> cigar>> mchr>> mpos>> isize>> seq>> qual;
//     EXPECT_STREQ(fragmentTag2.c_str()+1, tag.c_str());
//     EXPECT_STREQ("13", flag.c_str());
//     EXPECT_STREQ("1", S(chr));
//     EXPECT_STREQ("1", S(pos));
//     EXPECT_STREQ("0" , S(mapq));
//     EXPECT_STREQ("35M" , S(cigar));
//     EXPECT_STREQ("=" , S(mchr));
//     EXPECT_STREQ("1" , S(mpos));
//     EXPECT_STREQ("35" , S(isize));
//     EXPECT_STREQ(S(readFragment2), S(seq));
//     EXPECT_STREQ(S(dataQuality) , S(qual));
// }

// // TODO
// // Test color space, single end, 1 insertion
// // Test color space, single end, 1 deletion

#endif
