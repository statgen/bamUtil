#include <gtest/gtest.h>

#include "GenomeSequence.h"
#include "MatchedReadBase.h"
/*
 *  Explanatin of a position relationship graph
 *  Forward match:
 *  reference index:        0 1 2 3 4 5 6 7 8
 *  reference (BS):         A C T G A C G T A
 *  reference (CS):           1 2 1 2 1 3 1 3 
 *  read (CS):            A 0|1 2 1 2               (note: 1 2 1 2 is the part that actually aligned)
 *  read (BS):           (A)A|C T G A               (note: (A) is a primer. it should not be outputed)
 *  read index(CS):       0 1|2 3 4 5 
 *  read index(BS):         0|1 2 3 4
 *  genomeMatchPosition:     |^                     (note: ^ marked the genomeMatchPosition)          
 *  
 *  Backward match:
 *  reference index:        0 1 2 3 4 5 6 7 8
 *  reference (BS):         A C T G A C G T A
 *  reference (CS):           1 2 1 2 1 3 1 3 
 *  read (CS):                        1 3 1 3|0 T     (note: 2 1 3 1 is the part that actually aligned)
 *  read (BS - observed.)             T G C A|T(T)    (note: (T) is a primer, should not be outputed)
 *  read (BS - rev. comp.):           A C G T|A(A)    (note: reverse complement in base space) 
 *  read index (CS):                  5 4 3 2|1 0
 *  genomeMatchPosition:              ^               (note: ^ marked the genomeMatchPosition)          
 */

class MatchedReadBaseTest: public ::testing::Test
{
public:
    // const char* gs   = "ACTGACGTA";
    // const char* csgs = "512121313";
    GenomeSequence* gs;
    GenomeSequence* csgs;

    std::string cs_fwd_read;
    std::string cs_fwd_qual;
    std::string bs_fwd_read;
    std::string bs_fwd_qual; 
    std::string bs_fwd_read_correct;
    std::string bs_fwd_qual_correct;
    std::string cs_bwd_read;
    std::string cs_bwd_qual;
    std::string bs_bwd_read;
    std::string bs_bwd_qual; 
    std::string bs_bwd_read_correct;
    std::string bs_bwd_qual_correct;
    MatchedReadBase m;

protected:
    virtual void SetUp(){
        gs = new GenomeSequence;
        csgs = new GenomeSequence;
        gs->setReferenceName("short.fa");
        gs->open(false);
        csgs->setReferenceName("short.fa");
        csgs->open(true);

        cs_fwd_read = "A01212";
        cs_fwd_qual = "!12345";
        bs_fwd_read = "AACTGA";
        bs_fwd_qual = "!12345"; 

        bs_fwd_read_correct = bs_fwd_read;
        bs_fwd_qual_correct = bs_fwd_qual;

        cs_bwd_read = "T03131";
        cs_bwd_qual = "!12345";
        bs_bwd_read = "TTACGT";
        bs_bwd_qual = "!12345"; 

        bs_bwd_read_correct = bs_bwd_read;
        bs_bwd_qual_correct = bs_bwd_qual;
    }
    virtual void TearDown() {
        gs->close();
        delete gs;
        csgs->close();
        delete csgs;
    }
};

TEST_F(MatchedReadBaseTest, fixBaseRange_fwd_0base)
{
    uint32_t start = 2, end = 2;
    genomeIndex_t referencePos = 1;
    // not correct any base
    m.fixBaseRange(start, end,
                   bs_fwd_read, bs_fwd_qual,
                   cs_fwd_read, cs_fwd_qual,
                   gs, csgs, referencePos,
                   true);
    EXPECT_EQ(bs_fwd_read, bs_fwd_read_correct);
    EXPECT_EQ(bs_fwd_qual, bs_fwd_qual_correct);
}

TEST_F(MatchedReadBaseTest, fixBaseRange_fwd_1base)
{
    uint32_t start = 2, end = 2;
    genomeIndex_t referencePos = 1;
    // correct forward match, let us change 1 base
    cs_fwd_read[2]='3'; 
    bs_fwd_qual_correct[2] = '!' + 1;
    m.fixBaseRange(start, end,
                   bs_fwd_read, bs_fwd_qual,
                   cs_fwd_read, cs_fwd_qual,
                   gs, csgs, referencePos,
                   true);
    EXPECT_EQ(bs_fwd_read, bs_fwd_read_correct);
    EXPECT_EQ(bs_fwd_qual, bs_fwd_qual_correct);
}

TEST_F(MatchedReadBaseTest, fixBaseRange_fwd_2base)
{
    uint32_t start = 3, end = 4;
    genomeIndex_t referencePos = 1;
    // correct forward match, let us change 2 base
    // 2 base represents a SNP.
    // we change from 12 to 21
    cs_fwd_read[start]='1'; 
    cs_fwd_read[end]='2'; 
    bs_fwd_read_correct[start]='A'; // changed from 'T';
    m.fixBaseRange(start, end,
                   bs_fwd_read, bs_fwd_qual,
                   cs_fwd_read, cs_fwd_qual,
                   gs, csgs, referencePos,
                   true);
    EXPECT_EQ(bs_fwd_read, bs_fwd_read_correct);
    EXPECT_EQ(bs_fwd_qual, bs_fwd_qual_correct);

    // 2 base does not represent SNP.
    cs_fwd_read[start]='2'; 
    cs_fwd_read[end]='2'; 
    m.fixBaseRange(start, end,
                   bs_fwd_read, bs_fwd_qual,
                   cs_fwd_read, cs_fwd_qual,
                   gs, csgs, referencePos,
                   true); // will debug here
    // the following magic number 1 and 2 are hand picked.
    bs_fwd_qual_correct[start] = '!' + 1;
    bs_fwd_qual_correct[end] = '!' + 2;
    EXPECT_EQ(bs_fwd_read, bs_fwd_read_correct);
    EXPECT_EQ(bs_fwd_qual, bs_fwd_qual_correct);
}

TEST_F(MatchedReadBaseTest, fixBaseRange_fwd_toEnd)
{
    uint32_t start = 4, end = 5;
    genomeIndex_t referencePos = 1;
    // correct forward match, let us change 1 base
    cs_fwd_read[start]='0'; 
    cs_fwd_read[end]='0'; 
    bs_fwd_read_correct[start] = bs_fwd_read_correct[start-1];
    bs_fwd_read_correct[end] = bs_fwd_read_correct[end-1];

    m.fixBaseRange(start, end,
                   bs_fwd_read, bs_fwd_qual,
                   cs_fwd_read, cs_fwd_qual,
                   gs, csgs, referencePos,
                   true);
    EXPECT_EQ(bs_fwd_read, bs_fwd_read_correct);
    EXPECT_EQ(bs_fwd_qual, bs_fwd_qual_correct);
}


TEST_F(MatchedReadBaseTest, fixBaseRange_bwd_0base)
{
    uint32_t start = 2, end = 2;
    genomeIndex_t referencePos = 5;
    // not correct any base
    m.fixBaseRange(start, end,
                   bs_bwd_read, bs_bwd_qual,
                   cs_bwd_read, cs_bwd_qual,
                   gs, csgs, referencePos,
                   false);
    EXPECT_EQ(bs_bwd_read, bs_bwd_read_correct);
    EXPECT_EQ(bs_bwd_qual, bs_bwd_qual_correct);
}

TEST_F(MatchedReadBaseTest, fixBaseRange_bwd_1base)
{
    uint32_t start = 2, end = 2;
    genomeIndex_t referencePos = 5;
    // correct forward match, let us change 1 base
    cs_bwd_read[2]='0'; 
    bs_bwd_qual_correct[2] = '!' + 1;
    m.fixBaseRange(start, end,
                   bs_bwd_read, bs_bwd_qual,
                   cs_bwd_read, cs_bwd_qual,
                   gs, csgs, referencePos,
                   false);
    EXPECT_EQ(bs_bwd_read, bs_bwd_read_correct);
    EXPECT_EQ(bs_bwd_qual, bs_bwd_qual_correct);
}

TEST_F(MatchedReadBaseTest, fixBaseRange_bwd_2base)
{
    uint32_t start = 3, end = 4;
    genomeIndex_t referencePos = 5;
    // correct forward match, let us change 2 base
    // 2 base represents a SNP.
    // we change from 13 to 31
    cs_bwd_read[start]='3'; 
    cs_bwd_read[end]='1'; 
    bs_bwd_read_correct[start]='T'; // changed from 'A';
    m.fixBaseRange(start, end,
                   bs_bwd_read, bs_bwd_qual,
                   cs_bwd_read, cs_bwd_qual,
                   gs, csgs, referencePos,
                   false);
    EXPECT_EQ(bs_bwd_read, bs_bwd_read_correct);
    EXPECT_EQ(bs_bwd_qual, bs_bwd_qual_correct);

    // 2 base does not represent SNP.
    cs_bwd_read[start]='2'; 
    cs_bwd_read[end]='1'; 
    m.fixBaseRange(start, end,
                   bs_bwd_read, bs_bwd_qual,
                   cs_bwd_read, cs_bwd_qual,
                   gs, csgs, referencePos,
                   false); // will debug here
    // the following magic number 1 and 2 are hand picked.
    bs_bwd_qual_correct[start] = '!' + 1;
    bs_bwd_qual_correct[end] = '!' + 2;
    EXPECT_EQ(bs_bwd_read, bs_bwd_read_correct);
    EXPECT_EQ(bs_bwd_qual, bs_bwd_qual_correct);
}

TEST_F(MatchedReadBaseTest, fixBaseRange_bwd_toEnd)
{
    uint32_t start = 4, end = 5;
    genomeIndex_t referencePos = 5;
    // correct forward match, let us change 1 base
    cs_bwd_read[start]='0'; 
    cs_bwd_read[end]='0'; 
    bs_bwd_read_correct[start] = bs_bwd_read_correct[start-1];
    bs_bwd_read_correct[end] = bs_bwd_read_correct[end-1];

    m.fixBaseRange(start, end,
                   bs_bwd_read, bs_bwd_qual,
                   cs_bwd_read, cs_bwd_qual,
                   gs, csgs, referencePos,
                   false);
    EXPECT_EQ(bs_bwd_read, bs_bwd_read_correct);
    EXPECT_EQ(bs_bwd_qual, bs_bwd_qual_correct);
}

