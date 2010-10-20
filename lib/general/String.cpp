
#include <iostream>

#include "String.h"
#include <vector>

#ifdef OBSOLETE

std::vector<csg::string> *csg::string::split(char splitChar)
{
    std::vector<csg::string> *result = new std::vector<csg::string>;
    csg::string word;

    for (size_t i = 0; i<size(); i++)
    {
        if ((*this)[i]==splitChar)
        {
            result->push_back(word);
            word.clear();
        }
        else
            word.push_back((*this)[i]);
    }
    if (word.size()>0) result->push_back(word);
    return result;
}


#if defined(TEST)

int main(int argc, const char **argv)
{
    csg::string string("abcdef:abcdefghijk");

    std::vector<csg::string>    *result = string.split(':');

    for (int i=0; i<result->size(); i++)
    {
        std::cout << i << "\t" << (*result)[i] << std::endl;
    }
    delete result;  // suck

}
#endif

#endif
