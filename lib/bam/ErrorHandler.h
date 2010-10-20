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

#ifndef __ERROR_HANDLER_H__
#define __ERROR_HANDLER_H__

#include <iostream>

class ErrorHandler
{
public:

    // This specifies how this class should respond to errors.
    //   EXCEPTION - throw an exception for the error
    //   ABORT     - exit the program on the error
    //   RETURN    - just return failure on the error.
    enum HandlingType {EXCEPTION, ABORT, RETURN};

//     // Type of Error.
//     //    SUCCESS      : method completed successfully.
//     //    UNKNOWN      : unknown result (default value should never be used)
//     //    FAIL_IO      : method failed due to an I/O issue.
//     //    FAIL_MEM     : fail a memory allocation.
//     //    FAIL_ORDER   : method failed because it was called out of order,
//     //                   like trying to read a file without opening it for
//     //                   read or trying to read a record before the header.
//     //    FAIL_PARSE   : failed to parse a record/header - invalid format.
//     //    INVALID      : record is invalid other than for sorting.
//     //    NO_MORE_RECS : failed to read a record since there are no more to read
//     //                   either in the file or section if section based reading.
//     //    INVALID_SORT : record is invalid due to it not being sorted.
//     enum Type {SUCCESS = 0, UNKNOWN, FAIL_IO, FAIL_MEM, FAIL_ORDER, FAIL_PARSE,
//                INVALID, NO_MORE_RECS, INVALID_SORT};


    // Constructor
    ErrorHandler();
   
    // Destructor
    ~ErrorHandler();

    // Handle an error based on the error handling type.
    static void handleError(const char* message, 
                            HandlingType handlingType = EXCEPTION);
      
private:
};


#endif
