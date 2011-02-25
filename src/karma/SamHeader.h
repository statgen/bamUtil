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

#ifndef SAM_HEADER_H
#define SAM_HEADER_H

#include <string>
#include <iostream>
#include <map>

class SamHeader
{
private:
    std::map<std::string, std::map<std::string, std::string> > header;
public:
    std::map<std::string, std::string> &operator[](std::string &key)
    {
        return header[key];
    }
    std::map<std::string, std::string> &operator[](const char* key)
    {
        std::string k = key;
        return header[k];
    }
    void set(std::string &key1, std::string &key2, std::string &val);
    void set(const char *key1, const char *key2, const char *val)
    {
        std::string k1=key1, k2=key2, v=val;
        set(k1,k2,v);
    }
    void dump(std::ostream &);
    void clear() { header.clear(); };
    bool conformSpecification(); // check if the header conforms to SAM specification
    bool containTag(const std::string& key) { 
        if (header.find(key) != header.end()) return true;
        return false;
    };
    bool containTag(const char* key) { 
        std::string k = key;
        return containTag(k);
    };

};

#endif
