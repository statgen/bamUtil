#include <string>
#include <TrimSequence.h>

#include <gtest/gtest.h>

TEST(TrimSequenceTest, trimSequenceTest)
{
    std::string test;
    std::string::iterator result;

    test = "445566";
    result = trimSequence(test, '5', true);
    EXPECT_EQ(result - test.begin() , 2);

    test = "445554555";
    result = trimSequence(test, '5', true);
    EXPECT_EQ(result - test.begin(), 6);

    test = "4455545556";
    result = trimSequence(test, '5', true);
    EXPECT_EQ(result - test.begin(), 6);

    test = "44555455566";
    result = trimSequence(test, '5', true);
    EXPECT_EQ(result - test.begin(), 6);

    test = "665544";
    result = trimSequence(test, '5', false);
    EXPECT_EQ(test.end() - result , 2);

    test = "555455544";
    result = trimSequence(test, '5', false);
    EXPECT_EQ(test.end() - result, 6);

    test = "6555455544";
    result = trimSequence(test, '5', false);
    EXPECT_EQ(test.end() - result, 6);

    // Paul's test cases in TrimSequence.cpp 
    //
    // from the left:
    //
    test = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'A', true);
    EXPECT_TRUE(result == test.begin());

    test = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, '~', true);
    EXPECT_TRUE(result == test.end());

    test = "AAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'B', true);
    EXPECT_TRUE(result == (test.begin() + 5));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'B', true);
    EXPECT_TRUE(result == (test.begin() + 8));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'F', true);
    EXPECT_TRUE(result == (test.begin() + 12));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, '@', true);
    EXPECT_TRUE(result == (test.begin() + 0));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, '@', true);
    EXPECT_TRUE(result == (test.begin() + 0));

    test = "AAAFAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'F', true);
    EXPECT_TRUE(result == (test.begin() + 12));

    //
    // from the right:
    //
    test = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
    result = trimSequence(test, 'A', false);
    EXPECT_TRUE(result == test.end());

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
    result = trimSequence(test, '~', false);
    EXPECT_TRUE(result == test.begin());

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAA";
    result = trimSequence(test, 'B', false);
    EXPECT_TRUE(result == (test.end() - 5));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAAAA";
    result = trimSequence(test, 'B', false);
    EXPECT_TRUE(result == (test.end() - 7));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAAAAA";
    result = trimSequence(test, 'F', false);
    EXPECT_TRUE(result == (test.end() - 12));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAAAAA";
    result = trimSequence(test, '@', false);
    EXPECT_TRUE(result == (test.end() + 0));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAFAAA";
    result = trimSequence(test, 'F', false);
    EXPECT_TRUE(result == (test.end() - 12));
};
