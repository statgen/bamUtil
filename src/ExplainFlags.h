/*
 *  Copyright (C) 2013  Regents of the University of Michigan
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

#ifndef __EXPLAIN_FLAGS_H__
#define __EXPLAIN_FLAGS_H__

#include "BamExecutable.h"
#include "SamFile.h"

class ExplainFlags : public BamExecutable
{
public:
    ExplainFlags();
    static void explainFlagsDescription();
    void description();
    void usage();
    int execute(int argc, char **argv);
    virtual const char* getProgramName() {return("bam:explainFlags");}

private:
};

#endif
