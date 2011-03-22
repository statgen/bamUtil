/*
 * Copyright (c) 2009 Regents of the University of Michigan
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SamHeader.h"

void SamHeader::dump(std::ostream &f)
{
    std::map<std::string, std::map<std::string, std::string> >::iterator h;
    std::map<std::string, std::string>::iterator val;

    for (h = header.begin(); h != header.end(); h++)
    {
        f << "@" << h->first;
//        fprintf(f,"@%s", h->first);
        for (val=h->second.begin(); val!=h->second.end(); val++)
        {
            f << "\t" << val->first << ":" <<  val->second;
//            fprintf(f,"\t%s:%s", val->first.c_str(), val->second.c_str());
        }
        f << std::endl;
//        fprintf(f,"\n");
    }
}

void SamHeader::set(std::string &k1, std::string &k2, std::string &val)
{
    std::map<std::string, std::map<std::string, std::string> >::iterator h;

    h = header.find(k1);
    if (h==header.end())
    {
        std::map<std::string, std::string> newvals;
        header[k1] = newvals;
    }
    // gotta bomb here:
    header[k1][k2] = val;
}
