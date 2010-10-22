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

#ifndef __GLF_EXCEPTION_H__
#define __GLF_EXCEPTION_H__

#include <stdexcept> // stdexcept header file

#include "GlfStatus.h"

// GlfException objects should be thrown by functions that operate on 
// Glf files.
class GlfException : public std::exception
{
public:
    // constructor specifies default error message
    GlfException();
    GlfException(const std::string& what_arg);
    GlfException(GlfStatus::Status status, const std::string& errorMsg);
    GlfException(const GlfStatus& status);
    virtual ~GlfException() throw();

    virtual const char* what() const throw();

private:
    GlfStatus myStatus;
}; // end class GlfException


#endif
