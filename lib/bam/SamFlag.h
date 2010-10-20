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

#ifndef __SAM_FLAG_H__
#define __SAM_FLAG_H__

#include <stdint.h>


class SamFlag
{
public:
    /// Constants for parsing a flag.
    static const int16_t PAIRED              = 0x0001;
    static const int16_t PROPER_PAIR         = 0x0002;
    static const int16_t UNMAPPED            = 0x0004;
    static const int16_t MATE_UNMAPPED       = 0x0008;
    static const int16_t REVERSE             = 0x0010;
    static const int16_t MATE_REVERSED       = 0x0020;
    static const int16_t FIRST_READ          = 0x0040;
    static const int16_t SECOND_READ         = 0x0080;
    static const int16_t SECONDARY_ALIGNMENT = 0x0100;
    static const int16_t FAILED_QUALITY      = 0x0200;
    static const int16_t DUPLICATE           = 0x0400;

    static inline bool isMapped(uint16_t flag) {return(!(flag & UNMAPPED));}
    static inline bool isPaired(uint16_t flag) {return(flag & PAIRED);}
    static inline bool isReverse(uint16_t flag) {return(flag & REVERSE);}
    static inline bool isProperPair(uint16_t flag) {return(flag & PROPER_PAIR);}
    static inline bool isDuplicate(uint16_t flag) {return(flag & DUPLICATE);}
    static inline bool isQCFailure(uint16_t flag) {return(flag & FAILED_QUALITY);}

    static inline void setUnmapped(uint16_t& flag) { flag |= UNMAPPED;}

private:
    SamFlag();
};


#endif
