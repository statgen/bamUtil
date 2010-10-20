#ifndef _CSG_STRING_H
#define _CSG_STRING_H

#include <vector>

#pragma message "File String.h is obsolete... use STLUtilities.h instead"

#ifdef OBSOLETE
namespace csg
{

// provide a class derived from std::string to
// prevent excess rewriting of code.
//
class string : public std::string
{
public:
    string()
    {
        ;
    }
    string(const char *s) : std::string(s)
    {
        ;
    }

    // slowly split a string into a vector of strings
    std::vector<csg::string>    *split(char splitChar='\t');
};
};
#endif

#endif
