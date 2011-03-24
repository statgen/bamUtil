#include "ColorSpace.h"

TEST(MapperBaseTest, ColorSpaceRelatedFunction)
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

