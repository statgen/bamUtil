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

#ifndef __BAM_EXECUTABLE_H__
#define __BAM_EXECUTABLE_H__

#include "StringBasics.h"
#include "Parameters.h"

/// Base Class BAM Executable.
class BamExecutable
{
public: 
    static void bamVersion();
    static void bamExecutableDescription();
    BamExecutable();
    virtual ~BamExecutable();

    /// Print the 
    virtual void description();
    virtual void usage();
    virtual int execute(int argc, char**argv) = 0;
    
    virtual const char* getProgramName() {return("bam");}


protected:

private:
};

#endif
