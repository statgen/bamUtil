#include "GenomeSequence.h"

#include <algorithm>

#include "ColorSpace.h"

//////////////////////////////////////////////////////////////////////
// functor class for get the compliment base
//////////////////////////////////////////////////////////////////////
struct BaseComplementor
{
    char operator()(char& c)
    {
        return GenomeSequence::base2complement[(int)c];
    };
    char operator()(char c)
    {
        return GenomeSequence::base2complement[(int)c];
    };
};

// color space related variables and functions
///
///Convert 1 base and 1 color to 1 base (decode color space)
///
void BaseAndColorToBase(const char& base, const char& color, char& secondBase)
{
    int c = color - '0';
    BaseAndColorToBase(base, c, secondBase);
    return;
};

///
///Convert 1 base and 1 color to 1 base (decode color space)
///
void BaseAndColorToBase(const char& base, const int& color, char& secondBase)
{
    int hi = GenomeSequence::base2int[(int)(base)];
    int lo = color;
    if (hi>= 4 || lo >= 4 || lo < 0)
    {
        secondBase = 'N';
        return;
    }
    secondBase= GenomeSequence::int2base[(int)(fromCS2Base[hi<<2 | lo])];
    return;

};

///
///Convert 1 base and 1 base to the color of the transiation
///
void BaseAndBaseToColor(const char& base1, const char& base2, int& color)
{
    int hi = GenomeSequence::base2int[(int)(base1)];
    int lo = GenomeSequence::base2int[(int)(base2)];
    if (hi >= 4 || lo >= 4)
    {
        color = 5;
        return;
    }
    color = fromCS2Base[hi<<2 | lo];
    return;
};

///
///Convert 1 base and 1 base to the color of the transiation
///
void BaseAndBaseToColor(const char& base1, const char& base2, char& color)
{
    int c;
    BaseAndBaseToColor(base1, base2, c);
    color='0'+c;
    return;
};

/*
 * convert a color space read to base space
 * Warning: if the input read contains '5', the unknown base, then the output cannot be predicted.
 * @param keepPrimer true: the leading primer the colr space will be stored as the first base of the output, default value is [false]
 * @return translated read in base space
 */
std::string convertCSRead(const std::string& csread, bool keepPrimer)
{
    int cslength = csread.size();
    int bslength = (keepPrimer ? cslength : (cslength - 1));
    std::string bsread;
    bsread.resize(bslength);

    char lastBase = csread[0];
    bsread[0] = csread[0];
    for (int i= (keepPrimer ? 1: 0), j = 1;
         i < bslength; i++, j++)
    {
        BaseAndColorToBase(lastBase, csread[j], bsread[i]);
        lastBase = bsread[i];
    }
    return bsread;
};

std::string convertCSQuality(const std::string& csqual)
{
    std::string bsqual = csqual.substr(1);
    return bsqual;
};

/*
 * get reverse of a input read (in base space or in color space)
 * @param read input read
 * @return reverse of the input read
 */
std::string getReverse(const std::string& read)
{
    std::string ret = read;
    std::reverse(ret.begin(), ret.end());
    return ret;
};

/*
 * get compement of a input read (in base space or in color space)
 * @param read input read
 * @return complement of the input read
 */
std::string getComplement(const std::string& read)
{
    std::string ret = read;
    // here we can user functor+transform() to speed things up
    static BaseComplementor baseComplementor;
    std::transform(read.begin(), read.end(), ret.begin(), baseComplementor);
    return ret;
};

/*
 * get reverse complement of a input read (in base space or in color space)
 * @param read input read
 * @return reverse complement of the input read
 */
std::string getReverseComplement(const std::string& read)
{
    return getReverse(getComplement(read));
};

