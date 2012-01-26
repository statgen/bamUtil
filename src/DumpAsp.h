/*
 *  Copyright (C) 2012  Regents of the University of Michigan
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

//////////////////////////////////////////////////////////////////////////
// This file contains the processing for the executable option "dumpAsp"
// which prints the asp file to the screen in a readable format.

#ifndef __DUMP_ASP_H__
#define __DUMP_ASP_H__

#include "BamExecutable.h"
#include "AspRecord.h"

class DumpAsp : public BamExecutable
{
public:
    static void dumpAspDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);

private:
    void printRefOnly(AspRecord& record);
    void printDetailed(AspRecord& record);
};

#endif
