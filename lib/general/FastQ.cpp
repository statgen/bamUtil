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

#include "FastQ.h"

void FastQ::read(std::istream &stream)
{
    std::getline(stream, tag);
    std::getline(stream, data);
    std::getline(stream, tag2);
    std::getline(stream, quality);
    if (stream.eof())
    {
        clear();
    }
    else valid=true;
}

void FastQ::write(std::ostream &stream)
{
    if (valid)
    {
        stream << tag << std::endl;
        stream << data << std::endl;
        stream << tag2 << std::endl;
        stream << quality << std::endl;
    }
}


