/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
