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

#ifndef __SAM_TAGS_H__
#define __SAM_TAGS_H__

#include <stdint.h>


class SamTags
{
public:
    /// Constants for parsing a tags.
    static const char* BQ_TAG;
    static const char BQ_TAG_TYPE;
    static const char* ORIG_POS_TAG;
    static const char ORIG_POS_TAG_TYPE;
    static const char* ORIG_CIGAR_TAG;
    static const char ORIG_CIGAR_TAG_TYPE;
    static const char* ORIG_QUAL_TAG;
    static const char ORIG_QUAL_TAG_TYPE;

private:
    SamTags();
};


#endif
