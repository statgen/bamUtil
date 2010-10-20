#include "TrimSequence.h"

#if defined(TEST)

#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <string>

int main(int argc, const char **argv)
{
    std::string test;
    std::string::iterator result;

    //
    // from the left:
    //
    test = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'A', true);
    assert(result == test.begin());

    test = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, '~', true);
    assert(result == test.end());

    test = "AAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'B', true);
    assert(result == (test.begin() + 5));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'B', true);
    assert(result == (test.begin() + 8));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'F', true);
    assert(result == (test.begin() + 12));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, '@', true);
    assert(result == (test.begin() + 0));

    test = "AAAAAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, '@', true);
    assert(result == (test.begin() + 0));

    test = "AAAFAAAABCDEFGHIJKLMNOPQRSTUVWXYZ";
    result = trimSequence(test, 'F', true);
    assert(result == (test.begin() + 12));

    // trim left 12 bases, and untrimmed bases are 'FG' (turn bug into this test cass)
    test = "AAAFAAAABCDEFG";
    result = trimSequence(test, 'F', true);
    assert(result == (test.begin() + 12));

    //
    // from the right:
    //
    test = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
    result = trimSequence(test, 'A', false);
    assert(result == test.end());

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBA";
    result = trimSequence(test, '~', false);
    assert(result == test.begin());

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAA";
    result = trimSequence(test, 'B', false);
    assert(result == (test.end() - 5));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAAAA";
    result = trimSequence(test, 'B', false);
    assert(result == (test.end() - 7));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAAAAA";
    result = trimSequence(test, 'F', false);
    assert(result == (test.end() - 12));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAAAAA";
    result = trimSequence(test, '@', false);
    assert(result == (test.end() + 0));

    test = "ZYXWVUTSRQPONMLKJIHGFEDCBAAAAFAAA";
    result = trimSequence(test, 'F', false);
    assert(result == (test.end() - 12));

    test = "#################################";
    result = trimSequence(test, 'F', false);
    assert(result == (test.begin()));

#if 0
    // TODO: add explanation why this test case should trim 5 right most bases?
    test = ">BC@>28B==>=><?@=?>@8(>0309261/;6=@";
    result = trimSequence(test, '0', false);
    assert(result == (test.end())-5);
#endif

    exit(0);
}

#endif
