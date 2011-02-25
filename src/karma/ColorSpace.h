#ifndef _COLORSPACE_H_
#define _COLORSPACE_H_

#include <string>
// anything outside these values represents an invalid base
// base codes: 0-> A,    1-> C,     2-> G,      3-> T
// colorspace: 0-> blue, 1-> green, 2-> oragne, 3->red
//
// highest 2 bits: a base space nucleotide
// lowest 2 bits : a color space code
// Surprisingly, the following table is the same as const char fromBase2CS
static const char fromCS2Base[] =
{
    /* 0000 */ 0,   // A->A
    /* 0001 */ 1,   // A->C
    /* 0010 */ 2,   // A->G
    /* 0011 */ 3,   // A->T
    /* 0100 */ 1,   // C->C
    /* 0101 */ 0,   // C->A
    /* 0110 */ 3,   // C->T
    /* 0111 */ 2,   // C->G
    /* 1000 */ 2,   // G->G
    /* 1001 */ 3,   // G->T
    /* 1010 */ 0,   // G->A
    /* 1011 */ 1,   // G->C
    /* 1100 */ 3,   // T->T
    /* 1101 */ 2,   // T->G
    /* 1110 */ 1,   // T->C
    /* 1111 */ 0,   // T->A
};

// color space related variables and functions
///
///Convert 1 base and 1 color to 1 base (decode color space)
///
void BaseAndColorToBase(const char& base, const char& color, char& secondBase);

///
///Convert 1 base and 1 color to 1 base (decode color space)
///
void BaseAndColorToBase(const char& base, const int& color, char& secondBase);

///
///Convert 1 base and 1 base to the color of the transiation
///
void BaseAndBaseToColor(const char& base1, const char& base2, int& color);

///
///Convert 1 base and 1 base to the color of the transiation
///
void BaseAndBaseToColor(const char& base1, const char& base2, char& color);

/*
 * convert a color space read to base space
 * Warning: if the input read contains '5', the unknown base, then the output cannot be predicted.
 * @param keepPrimer true: the leading primer the colr space will be stored as the first base of the oupt
 * @return translated read in base space
 */
std::string convertCSRead(const std::string& csread, bool keepPrimer=false);

std::string convertCSQuality(const std::string& csqual);

/*
 * get reverse of a input read (in base space or in color space)
 * @param read input read
 * @return reverse of the input read
 */
std::string getReverse(const std::string& read);

/*
 * get compement of a input read (in base space or in color space)
 * @param read input read
 * @return complement of the input read
 */
std::string getComplement(const std::string& read);

/*
 * get reverse complement of a input read (in base space or in color space)
 * @param read input read
 * @return reverse complement of the input read
 */
std::string getReverseComplement(const std::string& read);



#endif /* _COLORSPACE_H_ */
