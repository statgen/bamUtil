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

#include "BaseUtilities.h"

#include <ctype.h>

bool BaseUtilities::isAmbiguous(char base)
{
    switch(base)
    {
        case 'N':
        case 'n':
        case '.':
            return(true);
            break;
        default:
            break;
    };

    // Not 'N', 'n', or '.', so return false.
    return(false);
}

bool BaseUtilities::areEqual(char base1, char base2)
{
    // If they are the same, return true.
    if(base1 == base2)
    {
        return(true);
    }
    // If one of the bases is '=', return true.
    if((base1 == '=') || (base2 == '='))
    {
        return(true);
    }

    // Check both in upercase.
    if(toupper(base1) == toupper(base2))
    {
        // same in upper case.
        return(true);
    }

    // The bases are different.
    return(false);
}


// Get phred base quality from the specified ascii quality.
uint8_t BaseUtilities::getPhredBaseQuality(char charQuality)
{
    return(charQuality - 33);
}
