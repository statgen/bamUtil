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

#ifndef _FASTQ_H
#define _FASTQ_H
#include <iostream>
#include <string>
#include <math.h>

class FastQ
{
public:
    bool valid;
    FastQ()
    {
        clear();
    }
    void clear()
    {
        valid=false;
    }
    std::string  tag;
    std::string  data;
    std::string  tag2;
    std::string  quality;
    void read(std::istream &);
    void write(std::ostream &);
};

inline std::istream &operator >> (std::istream &stream, FastQ &f)
{
    f.read(stream);
    return stream;
}

inline std::ostream &operator << (std::ostream &stream, FastQ &f)
{
    f.write(stream);
    return stream;
}

inline double phred2double(const char& c)
{
    double errorProb=pow(10,-double((c-'!'))/10);
    return errorProb;
}

#endif
