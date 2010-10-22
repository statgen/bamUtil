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

#include <iostream>
#include <string.h>
#include <stdlib.h>

#include "DumpGlf.h"

void Usage()
{
    std::cerr << "Usage: " << std::endl;
    std::cerr << std::endl;
    dumpGlfDescription();
}


int main(int argc, char ** argv)
{
    // Verify at least one arg.
    if(argc < 2)
    {
        // Not enough args...
        std::cerr << "Error: Not enough arguments\n";
        Usage();
        exit(-1);
    }

    if(strcmp(argv[1], "dumpGlf") == 0)
    {
        if(argc != 3)
        {
            dumpGlfUsage();
            exit(-1);
        }
        // Dump the glf file.
        return(dumpGlf(argv[2]));
    }
    else
    {
        std::cerr << "Error: Unknown Option\n";
        Usage();
        exit(-1);
    }
}
