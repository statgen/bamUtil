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

#include "WindowsHelper.h"
#ifdef    __WIN32__
#ifndef   __GNUC__
#include <dir.h>

void WildCardArguments(int & argc, char ** & argv)
{
    if (argc < 2) return;

    int  count = 0;
    for (int i = 1; i < argc; i++)
    {
        struct ffblk blk;

        int done = findfirst(argv[i], &blk, 0);
        while (!done)
        {
            done = findnext(&blk);
            count++;
        }
    }

    char ** new_argv = new char * [count + 1];
    int     new_argc = 1;

    new_argv[0] = argv[0];
    for (int i = 1; i < argc; i++)
    {
        struct ffblk blk;

        int done = findfirst(argv[i], &blk, 0);
        while (!done && new_argc <= count)
        {
            new_argv[new_argc++] = strdup(blk.ff_name);
            done = findnext(&blk);
        }
    }

    argc = new_argc;
    argv = new_argv;
}

#endif
#endif

