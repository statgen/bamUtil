

// put TEST below here, so that makedepend will see the .h, so that we
// can get a clean dependency for SmithWaterman.o, so that we can at least
// compile the header when we change it.


#include <getopt.h>
#include "Generic.h"
#include "CigarRollerTest.h"

//
// Test the obvious cases.
// Add non-obvious ones as bugs come up.
//
void CigarRollerTest::test(void)
{
    // Create the CigarRoller.
    CigarRoller cigar;

    int failures = 0, testNum = 0;

    //   const char *str;
    String str;
    std::string result;
    
    cigar.getCigarString(str);
    result = str.c_str();
    //    result = str = cigar.getString(); delete str;

    check(failures, ++testNum, "constructor", result, "");    // test empty case

    CigarRoller::CigarOperator op;

    op.operation = CigarRoller::match;
    op.count = 10;

    cigar += op;

    op.count=5;

    cigar += op;

    op.count=5; op.operation = CigarRoller::mismatch;   // test that match/mismatch get combined

    cigar += op;

    op.count=5; op.operation = CigarRoller::insert;

    cigar += op;

    op.count=5; op.operation = CigarRoller::insert;

    cigar += op;

    op.count=5; op.operation = CigarRoller::del;

    cigar += op;

    op.count=5; op.operation = CigarRoller::mismatch;

    cigar += op;

    op.count=5; op.operation = CigarRoller::match;

    cigar += op;

    op.count=5; op.operation = CigarRoller::skip;

    cigar += op;

    op.count=5; op.operation = CigarRoller::match;

    cigar += op;

    op.count=2; op.operation = CigarRoller::pad;

    cigar += op;

    op.count=3; op.operation = CigarRoller::match;

    cigar += op;

    cigar.getCigarString(str);
    result = str.c_str();
    //    result = str = cigar.getString();  delete str;

    check(failures, ++testNum, "match combining", "20M10I5D10M5N5M2P3M", result);
    check(failures, ++testNum, "length check", 8, cigar.size());

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Test getRefOffset, getQueryIndex, getRefPosition, & getQueryIndex(that takes ref position)

    // Validate the reference offsets when passing in a query index,
    // and the query offsets when passing in a query index.
    // Spot check the offsets out of order - to verify order doesn't matter.
    check(failures, ++testNum, "getRefOffset(20)", -1, cigar.getRefOffset(20));
    check(failures, ++testNum, "getRefOffset(30)", 25, cigar.getRefOffset(30));
    check(failures, ++testNum, "getRefOffset(46)", 46, cigar.getRefOffset(46));
    check(failures, ++testNum, "getRefOffset(0)", 0, cigar.getRefOffset(0));
    check(failures, ++testNum, "getRefPosition(20, 5)", -1, cigar.getRefPosition(20, 5));
    check(failures, ++testNum, "getRefPosition(30, 5)", 30, cigar.getRefPosition(30, 5));
    check(failures, ++testNum, "getRefPosition(46, 5)", 51, cigar.getRefPosition(46, 5));
    check(failures, ++testNum, "getRefPosition(0, 5)", 5, cigar.getRefPosition(0, 5));
    check(failures, ++testNum, "getQueryIndex(20)", CigarRoller::INDEX_NA, cigar.getQueryIndex(20));
    check(failures, ++testNum, "getQueryIndex(30)", 35, cigar.getQueryIndex(30));
    check(failures, ++testNum, "getQueryIndex(35)", Cigar::INDEX_NA, cigar.getQueryIndex(35));
    check(failures, ++testNum, "getQueryIndex(46)", 46, cigar.getQueryIndex(46));
    check(failures, ++testNum, "getQueryIndex(0)", 0, cigar.getQueryIndex(0));
    check(failures, ++testNum, "getQueryIndex(25, 5)", -1, cigar.getQueryIndex(20));
    check(failures, ++testNum, "getQueryIndex(35, 5)", 35, cigar.getQueryIndex(30));
    check(failures, ++testNum, "getQueryIndex(40, 5)", -1, cigar.getQueryIndex(35));
    check(failures, ++testNum, "getQueryIndex(51, 5)", 46, cigar.getQueryIndex(46));
    check(failures, ++testNum, "getQueryIndex(5, 5)", 0, cigar.getQueryIndex(0));

    int i = 0;
    int queryIndex = 0;
    int refOffset = 0;
    // Validate the 20 matches.
    for(i = 0; i < 20; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }
    // Validate the 10 Inserts - exist in query, but not in the reference.
    for(i = 0; i < 10; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", -1, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", -1, cigar.getRefPosition(queryIndex, 5));
        ++queryIndex;
    }
    // Validate the 5 Deletions - exist in the reference, but not the query.
    for(i = 0; i < 5; i++)
    {
        check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));
        refOffset++;
    }
    // Validate the 10 matches.
    for(i = 0; i < 10; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }
    // Validate the 5 Skips - exist in the reference, but not the query.
    for(i = 0; i < 5; i++)
    {
        check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));
        refOffset++;
    }
    // Validate the 5 matches.
    for(i = 0; i < 5; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }
    // Nothing to validate for the 2 pads since they do not show up in either the reference or the query.
    // Validate the 3 matches.
    for(i = 0; i < 3; i++)
    {
        check(failures, ++testNum, "getRefOffset(queryIndex)", refOffset, cigar.getRefOffset(queryIndex));
        check(failures, ++testNum, "getQueryIndex(refOffset)", queryIndex, cigar.getQueryIndex(refOffset));
        check(failures, ++testNum, "getRefPosition(queryIndex, 5)", refOffset + 5, cigar.getRefPosition(queryIndex, 5));
        check(failures, ++testNum, "getQueryIndex(refPosition, 5)", queryIndex, cigar.getQueryIndex(refOffset + 5, 5));
        ++queryIndex;
        ++refOffset;
    }

    // Validate that if you go beyond the end, -1 is returned.
    check(failures, ++testNum, "getRefOffset(queryIndex)", -1, cigar.getRefOffset(queryIndex));
    check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
    check(failures, ++testNum, "getRefPosition(queryIndex, 5)", -1, cigar.getRefPosition(queryIndex, 5));
    check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));
    ++queryIndex;
    ++refOffset;
    check(failures, ++testNum, "getRefOffset(queryIndex)", -1, cigar.getRefOffset(queryIndex));
    check(failures, ++testNum, "getQueryIndex(refOffset)", -1, cigar.getQueryIndex(refOffset));
    check(failures, ++testNum, "getRefPosition(queryIndex, 5)", -1, cigar.getRefPosition(queryIndex, 5));
    check(failures, ++testNum, "getQueryIndex(refPosition, 5)", -1, cigar.getQueryIndex(refOffset + 5, 5));

    ////////////////////////////////////////////////////////////////////////
    // Test getNumOverlaps
    
    // When query starts at position 5:
    // Overlaps are at 5-24, 30-39, 45-49, 50-52

    // Test with region [1-5) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(1, 5, 5)", (uint32_t)0, cigar.getNumOverlaps(1, 5, 5));

    // Test with region [53-101) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(53, 101, 5)", (uint32_t)0, cigar.getNumOverlaps(53, 101, 5));

    // Test with region [53-10) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(53, 10, 5)", (uint32_t)0, cigar.getNumOverlaps(53, 10, 5));

    // Test with region [35-10) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(35, 10, 5)", (uint32_t)0, cigar.getNumOverlaps(35, 10, 5));

    // Test with region [35-1) where query starts at position 5 = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(35, 1, 5)", (uint32_t)0, cigar.getNumOverlaps(35, 1, 5));

    // Test with region [1-6) where query starts at position 5 - 1 overlap = pos 5.
    check(failures, ++testNum, "getNumOverlaps(1, 6, 5)", (uint32_t)1, cigar.getNumOverlaps(1, 6, 5));

    // Test with region [25-30) where query has only DELETIONS from the reference = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(25, 30, 5)", (uint32_t)0, cigar.getNumOverlaps(25, 30, 5));

    // Test with region [24-30) where query has only a match at position 24 = 1 overlap.
    check(failures, ++testNum, "getNumOverlaps(24, 30, 5)", (uint32_t)1, cigar.getNumOverlaps(24, 30, 5));

    // Test with region [25-31) where query has only a match at position 30 = 1 overlap.
    check(failures, ++testNum, "getNumOverlaps(25, 31, 5)", (uint32_t)1, cigar.getNumOverlaps(25, 31, 5));

    // Test with region [1-31), match pos 5-24 & 30 = 21 overlaps
    check(failures, ++testNum, "getNumOverlaps(1, 31, 5)", (uint32_t)21, cigar.getNumOverlaps(1, 31, 5));

    // Test a region that covers the entire range [1-101), match pos 5-24, 30-39, 45-49, & 50-52 = 38 overlaps
    check(failures, ++testNum, "getNumOverlaps(1, 101, 5)", (uint32_t)38, cigar.getNumOverlaps(1, 101, 5));

    // Test a region that covers the entire range [-1--1), (whole region) match pos 5-24, 30-39, 45-49, & 50-52 = 38 overlaps
    check(failures, ++testNum, "getNumOverlaps(-1, -1, 5)", (uint32_t)38, cigar.getNumOverlaps(-1, -1, 5));

    // Test a region that covers the entire range [6-52), match pos 6-24, 30-39, 45-49, & 50-51 = 36 overlaps
    check(failures, ++testNum, "getNumOverlaps(6, 52, 5)", (uint32_t)36, cigar.getNumOverlaps(6, 52, 5));

    // Test with region [40-45) where query has only SKIPS from the reference = 0 overlaps.
    check(failures, ++testNum, "getNumOverlaps(40, 45, 5)", (uint32_t)0, cigar.getNumOverlaps(40, 45, 5));

    // Test a region that covers the range [-1-10), (whole region) match pos 5-9 = 5 overlaps
    check(failures, ++testNum, "getNumOverlaps(-1, 10, 5)", (uint32_t)5, cigar.getNumOverlaps(-1, 10, 5));

    // Test a region that covers the range [50--1), (whole region) match pos 50-52 = 3 overlaps
    check(failures, ++testNum, "getNumOverlaps(50, -1, 5)", (uint32_t)3, cigar.getNumOverlaps(50, -1, 5));

    ////////////////////////////////////////////////////////////////////////////
    // Test a new CIGAR.
    cigar.Set("4M10N4M3I2M4D3M");
    String expectedResult = "4M10N4M3I2M4D3M";
    String cigarString = "HI";
    cigar.getCigarString(cigarString);
    check(failures, ++testNum, "getCigarString", expectedResult, cigarString);

    ////////////////////////
    // Test getNumOverlaps.
    check(failures, ++testNum, "getNumOverlaps(5,32,5)", (uint32_t)13, cigar.getNumOverlaps(5,32,5));
    check(failures, ++testNum, "getNumOverlaps(5,31,5)", (uint32_t)12, cigar.getNumOverlaps(5,31,5));
    check(failures, ++testNum, "getNumOverlaps(0,100,5)", (uint32_t)13, cigar.getNumOverlaps(0,100,5));
    check(failures, ++testNum, "getNumOverlaps(-1, -1,5)", (uint32_t)13, cigar.getNumOverlaps(-1, -1,5));
    check(failures, ++testNum, "getNumOverlaps(-1,10,5)", (uint32_t)4, cigar.getNumOverlaps(-1,10,5));
    check(failures, ++testNum, "getNumOverlaps(10,-1,5)", (uint32_t)9, cigar.getNumOverlaps(10,-1,5));
    check(failures, ++testNum, "getNumOverlaps(9,19,5)", (uint32_t)0, cigar.getNumOverlaps(9,19,5));
    check(failures, ++testNum, "getNumOverlaps(9,20,5)", (uint32_t)1, cigar.getNumOverlaps(9,20,5));
    check(failures, ++testNum, "getNumOverlaps(9,6,5)", (uint32_t)0, cigar.getNumOverlaps(9,6,5));
    check(failures, ++testNum, "getNumOverlaps(0,5,5)", (uint32_t)0, cigar.getNumOverlaps(0,5,5));
    check(failures, ++testNum, "getNumOverlaps(32,40,5)", (uint32_t)0, cigar.getNumOverlaps(32,40,5));
    check(failures, ++testNum, "getNumOverlaps(0,5,1)", (uint32_t)4, cigar.getNumOverlaps(0,5,1));
    check(failures, ++testNum, "getNumOverlaps(32,40,32)", (uint32_t)4, cigar.getNumOverlaps(32,40,32));

    // Get Query Index for reference offset 0 - 27
    // 4M
    check(failures, ++testNum, "getQueryIndex(0)", 0, cigar.getQueryIndex(0));
    check(failures, ++testNum, "getQueryIndex(1)", 1, cigar.getQueryIndex(1));
    check(failures, ++testNum, "getQueryIndex(2)", 2, cigar.getQueryIndex(2));
    check(failures, ++testNum, "getQueryIndex(3)", 3, cigar.getQueryIndex(3));
    // 10N
    check(failures, ++testNum, "getQueryIndex(4)", -1, cigar.getQueryIndex(4));
    check(failures, ++testNum, "getQueryIndex(5)", -1, cigar.getQueryIndex(5));
    check(failures, ++testNum, "getQueryIndex(6)", -1, cigar.getQueryIndex(6));
    check(failures, ++testNum, "getQueryIndex(7)", -1, cigar.getQueryIndex(7));
    check(failures, ++testNum, "getQueryIndex(8)", -1, cigar.getQueryIndex(8));
    check(failures, ++testNum, "getQueryIndex(9)", -1, cigar.getQueryIndex(9));
    check(failures, ++testNum, "getQueryIndex(10)", -1, cigar.getQueryIndex(10));
    check(failures, ++testNum, "getQueryIndex(11)", -1, cigar.getQueryIndex(11));
    check(failures, ++testNum, "getQueryIndex(12)", -1, cigar.getQueryIndex(12));
    check(failures, ++testNum, "getQueryIndex(13)", -1, cigar.getQueryIndex(13));
    // 4M
    check(failures, ++testNum, "getQueryIndex(14)", 4, cigar.getQueryIndex(14));
    check(failures, ++testNum, "getQueryIndex(15)", 5, cigar.getQueryIndex(15));
    check(failures, ++testNum, "getQueryIndex(16)", 6, cigar.getQueryIndex(16));
    check(failures, ++testNum, "getQueryIndex(17)", 7, cigar.getQueryIndex(17));
    // 3I - nothing to check - not in reference - covers query indices  8-10
    // 2M
    check(failures, ++testNum, "getQueryIndex(18)", 11, cigar.getQueryIndex(18));
    check(failures, ++testNum, "getQueryIndex(19)", 12, cigar.getQueryIndex(19));
    // 4D
    check(failures, ++testNum, "getQueryIndex(20)", -1, cigar.getQueryIndex(20));
    check(failures, ++testNum, "getQueryIndex(21)", -1, cigar.getQueryIndex(21));
    check(failures, ++testNum, "getQueryIndex(22)", -1, cigar.getQueryIndex(22));
    check(failures, ++testNum, "getQueryIndex(23)", -1, cigar.getQueryIndex(23));
    // 3M
    check(failures, ++testNum, "getQueryIndex(24)", 13, cigar.getQueryIndex(24));
    check(failures, ++testNum, "getQueryIndex(25)", 14, cigar.getQueryIndex(25));
    check(failures, ++testNum, "getQueryIndex(26)", 15, cigar.getQueryIndex(26));

    //  Get Query Index for reference positions 0-33
    // N/A
    check(failures, ++testNum, "getQueryIndex(0, 5)", -1, cigar.getQueryIndex(0, 5));
    check(failures, ++testNum, "getQueryIndex(1, 5)", -1, cigar.getQueryIndex(1, 5));
    check(failures, ++testNum, "getQueryIndex(2, 5)", -1, cigar.getQueryIndex(2, 5));
    check(failures, ++testNum, "getQueryIndex(3, 5)", -1, cigar.getQueryIndex(3, 5));
    check(failures, ++testNum, "getQueryIndex(4, 5)", -1, cigar.getQueryIndex(4, 5));
    // 4M
    check(failures, ++testNum, "getQueryIndex(5, 5)", 0, cigar.getQueryIndex(5, 5));
    check(failures, ++testNum, "getQueryIndex(6, 5)", 1, cigar.getQueryIndex(6, 5));
    check(failures, ++testNum, "getQueryIndex(7, 5)", 2, cigar.getQueryIndex(7, 5));
    check(failures, ++testNum, "getQueryIndex(8, 5)", 3, cigar.getQueryIndex(8, 5));
    // 10N
    check(failures, ++testNum, "getQueryIndex(9, 5)", -1, cigar.getQueryIndex(9, 5));
    check(failures, ++testNum, "getQueryIndex(10, 5)", -1, cigar.getQueryIndex(10, 5));
    check(failures, ++testNum, "getQueryIndex(11, 5)", -1, cigar.getQueryIndex(11, 5));
    check(failures, ++testNum, "getQueryIndex(12, 5)", -1, cigar.getQueryIndex(12, 5));
    check(failures, ++testNum, "getQueryIndex(13, 5)", -1, cigar.getQueryIndex(13, 5));
    check(failures, ++testNum, "getQueryIndex(14, 5)", -1, cigar.getQueryIndex(14, 5));
    check(failures, ++testNum, "getQueryIndex(15, 5)", -1, cigar.getQueryIndex(15, 5));
    check(failures, ++testNum, "getQueryIndex(16, 5)", -1, cigar.getQueryIndex(16, 5));
    check(failures, ++testNum, "getQueryIndex(17, 5)", -1, cigar.getQueryIndex(17, 5));
    check(failures, ++testNum, "getQueryIndex(18, 5)", -1, cigar.getQueryIndex(18, 5));
    // 4M
    check(failures, ++testNum, "getQueryIndex(19, 5)", 4, cigar.getQueryIndex(19, 5));
    check(failures, ++testNum, "getQueryIndex(20, 5)", 5, cigar.getQueryIndex(20, 5));
    check(failures, ++testNum, "getQueryIndex(21, 5)", 6, cigar.getQueryIndex(21, 5));
    check(failures, ++testNum, "getQueryIndex(22, 5)", 7, cigar.getQueryIndex(22, 5));
    // 3I - nothing to check - not in reference - covers query indices  8-10
    // 2M
    check(failures, ++testNum, "getQueryIndex(23, 5)", 11, cigar.getQueryIndex(23, 5));
    check(failures, ++testNum, "getQueryIndex(24, 5)", 12, cigar.getQueryIndex(24, 5));
    // 4D
    check(failures, ++testNum, "getQueryIndex(25, 5)", -1, cigar.getQueryIndex(25, 5));
    check(failures, ++testNum, "getQueryIndex(26, 5)", -1, cigar.getQueryIndex(26, 5));
    check(failures, ++testNum, "getQueryIndex(27, 5)", -1, cigar.getQueryIndex(27, 5));
    check(failures, ++testNum, "getQueryIndex(28, 5)", -1, cigar.getQueryIndex(28, 5));
    // 3M
    check(failures, ++testNum, "getQueryIndex(29, 5)", 13, cigar.getQueryIndex(29, 5));
    check(failures, ++testNum, "getQueryIndex(30, 5)", 14, cigar.getQueryIndex(30, 5));
    check(failures, ++testNum, "getQueryIndex(31, 5)", 15, cigar.getQueryIndex(31, 5));

    // Get reference offset for query index 0 - 17
    // 4M
    check(failures, ++testNum, "getRefOffset(0)", 0, cigar.getRefOffset(0));
    check(failures, ++testNum, "getRefOffset(1)", 1, cigar.getRefOffset(1));
    check(failures, ++testNum, "getRefOffset(2)", 2, cigar.getRefOffset(2));
    check(failures, ++testNum, "getRefOffset(3)", 3, cigar.getRefOffset(3));
    // 10N - nothing to check - not in query - covers ref offsets 4-13
    // 4M
    check(failures, ++testNum, "getRefOffset(4)", 14, cigar.getRefOffset(4));
    check(failures, ++testNum, "getRefOffset(5)", 15, cigar.getRefOffset(5));
    check(failures, ++testNum, "getRefOffset(6)", 16, cigar.getRefOffset(6));
    check(failures, ++testNum, "getRefOffset(7)", 17, cigar.getRefOffset(7));
    // 3I
    check(failures, ++testNum, "getRefOffset(8)", -1, cigar.getRefOffset(8));
    check(failures, ++testNum, "getRefOffset(9)", -1, cigar.getRefOffset(9));
    check(failures, ++testNum, "getRefOffset(10)", -1, cigar.getRefOffset(10));
    // 2M
    check(failures, ++testNum, "getRefOffset(11)", 18, cigar.getRefOffset(11));
    check(failures, ++testNum, "getRefOffset(12)", 19, cigar.getRefOffset(12));
    // 4D - nothing to check - not in query - covers ref offsets 20-23
    // 3M
    check(failures, ++testNum, "getRefOffset(13)", 24, cigar.getRefOffset(13));
    check(failures, ++testNum, "getRefOffset(14)", 25, cigar.getRefOffset(14));
    check(failures, ++testNum, "getRefOffset(15)", 26, cigar.getRefOffset(15));


    // Get reference position for query index 0 - 17
    // 4M
    check(failures, ++testNum, "getRefPosition(0, 5)", 5, cigar.getRefPosition(0, 5));
    check(failures, ++testNum, "getRefPosition(1, 5)", 6, cigar.getRefPosition(1, 5));
    check(failures, ++testNum, "getRefPosition(2, 5)", 7, cigar.getRefPosition(2, 5));
    check(failures, ++testNum, "getRefPosition(3, 5)", 8, cigar.getRefPosition(3, 5));
    // 10N - nothing to check - not in query - covers ref offsets 4-13
    // 4M
    check(failures, ++testNum, "getRefPosition(4, 5)", 19, cigar.getRefPosition(4, 5));
    check(failures, ++testNum, "getRefPosition(5, 5)", 20, cigar.getRefPosition(5, 5));
    check(failures, ++testNum, "getRefPosition(6, 5)", 21, cigar.getRefPosition(6, 5));
    check(failures, ++testNum, "getRefPosition(7, 5)", 22, cigar.getRefPosition(7, 5));
    // 3I
    check(failures, ++testNum, "getRefPosition(8, 5)", -1, cigar.getRefPosition(8, 5));
    check(failures, ++testNum, "getRefPosition(9, 5)", -1, cigar.getRefPosition(9, 5));
    check(failures, ++testNum, "getRefPosition(10, 5)", -1, cigar.getRefPosition(10, 5));
    // 2M
    check(failures, ++testNum, "getRefPosition(11, 5)", 23, cigar.getRefPosition(11, 5));
    check(failures, ++testNum, "getRefPosition(12, 5)", 24, cigar.getRefPosition(12, 5));
    // 4D - nothing to check - not in query - covers ref pos 25-28
    // 3M
    check(failures, ++testNum, "getRefPosition(13, 5)", 29, cigar.getRefPosition(13, 5));
    check(failures, ++testNum, "getRefPosition(14, 5)", 30, cigar.getRefPosition(14, 5));
    check(failures, ++testNum, "getRefPosition(15, 5)", 31, cigar.getRefPosition(15, 5));



    ////////////////////////////////////////////////////////////////////////////
    // Test a new CIGAR set by buffer.
    // 2S 3M 1I 2M 1D 1M 2P 1M 3N 1M 3H
    uint32_t cigarBuf[] = {0x24,  // 2S = 2 << 4 | 4
                           0x30,  // 3M = 3 << 4 | 0
                           0x11,  // 1I = 1 << 4 | 1
                           0x20,  // 2M = 2 << 4 | 0
                           0x12,  // 1D = 1 << 4 | 2
                           0x10,  // 1M = 1 << 4 | 0
                           0x26,  // 2P = 2 << 4 | 6
                           0x10,  // 1m = 1 << 4 | 0
                           0x33,  // 3N = 3 << 4 | 3
                           0x10,  // 1M = 1 << 4 | 0
                           0x35}; // 3H = 3 << 4 | 5
    cigar.Set(cigarBuf, 11);
    cigarString = "HI";
    cigar.getCigarString(cigarString);
    expectedResult = "2S3M1I2M1D1M2P1M3N1M3H";
    check(failures, ++testNum, "getCigarString", expectedResult, cigarString);
    check(failures, ++testNum, "getNumEndClips", 3, cigar.getNumEndClips());
    check(failures, ++testNum, "getNumBeginClips", 2, cigar.getNumBeginClips());

    std::cout << "\nCigarRoller PASS: " << testNum - failures << "  FAIL: " << failures << std::endl;
}


int main(int argc, const char **argv)
{
    CigarRollerTest roller;
    
    bool showAllCasesFlag = false;
    int opt;

    while(( opt = getopt(argc, (char **) argv, "v")) != -1) {
        switch(opt) {
            case 'v':
                showAllCasesFlag = true;
                break;
            default:
                std::cerr << "usage: testSW [-v]" << std::endl;
                exit(1);
        }
    }


    //
    // do cigar roller tests first
    //
    roller.test();

    // CIGAR explanation - for backward SW runs, the corresponding
    // CIGAR string is generated from the back of the string to the
    // front.  Recall that the soft clipping is only done at the
    // "end" of the string, taking direction into account.

    // Comment out this result since it doesn't refelct the results of test.
    //    cout << endl << "Total Errors found: " << errors << endl;
}
