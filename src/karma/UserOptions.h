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

 #ifndef _USER_OPTIONS
 #define _USER_OPTIONS

 #ifndef __STDC_LIMIT_MACROS
 #define __STDC_LIMIT_MACROS
 #endif

 #include <stdint.h>

 #include <iostream>
 #include <stdint.h>
 #include <cstring>
 #include <cstdlib>
 #include <vector>
 #include <string>


 //
 // base class of a hierarchy of command argument processing.
 //
 // The goal is to minimize redundant argument checking and
 // help message printing, as well as to ensure consistent
 // argument processing among different modes.
 //
 // baseArguments is the base framework plus a few common
 // underlying flags.
 //
 // Derive other classes for other modes and add arguments
 // specific to the mode involved.
 //
class BaseArguments
{
  protected:
    std::string optionString;
    virtual void getoptH();
    int argc;
    const char **argv;
    int opt;
  public:
    bool debug;
    int seed;
  BaseArguments() : debug(false), seed(12345)
    {
        optionString += "ds:";
    }
    virtual ~BaseArguments() {}
    void setArgs(int argc, const char **argv)
    {
        this->argc = argc;
        this->argv = argv;
    }
    virtual void print(std::ostream &out)
    {
        out << "debug " << (debug ? "on" : "off") << " (default off)" << std::endl;
        out << "seed " << seed << " (default 12345)" << std::endl;
    }
    virtual void usage();  // pretty printing of values/defaults/etc
    void getopt();
};

class MapArguments : public BaseArguments
{
 public:
    int      mapInColorSpace;
    int      insertSize;
    uint64_t maxBases;
    uint64_t maxReads;
    int      occurrenceCutoff;
    int      qualityTrim;
    bool     showReferenceBases;// show all reference bases in the sam file
    bool     quietMode;           // quiet mode
    int      wordSize;
    int      numThread;
    std::string                 outputFilename;
    std::vector<std::string>    references;
    std::vector<std::string>    SAMHeaderOptions;
 MapArguments() :
    mapInColorSpace(false),
        insertSize(250),
        maxBases(0),
        maxReads(0),
        occurrenceCutoff(5000),
        qualityTrim(0),
        showReferenceBases(false),
        quietMode(false),
        wordSize(15),
        numThread(1)    
        {
            optionString += "a:B:cEH:m:o:O:q:r:R:t:Qw:";
        }
    ~MapArguments()
    {
        ;
    }
    void print(std::ostream &out)
    {
        out << "insertSize " << insertSize << " (default 250)" << std::endl;
        out << "SAMHeaderOptions:";
        for (std::vector<std::string>::iterator it=SAMHeaderOptions.begin(); it<SAMHeaderOptions.end(); it++)
        {
            out << " " << *it;
        }
        out << std::endl;
        out << "outputFilename " << outputFilename << " (no default)" << std::endl;
        out << "maxBases " << maxBases << " (default 0)" << std::endl;
        out << "maxReads " << maxReads << " (default 0)" << std::endl;
        out << "occurrenceCutoff " << occurrenceCutoff << " (default 5000)" << std::endl;
        out << "numThread " << numThread << " (default 1)" << std::endl;
        out << "quiet mode " << quietMode << " (default off)" << std::endl;;
        out << "show reference bases " << showReferenceBases << " (default off)" << std::endl;;
        out << "reference(s) (no default):";
        for (std::vector<std::string>::iterator it=references.begin(); it<references.end(); it++)
        {
            out << " " << *it;
        }
        out << std::endl;

        BaseArguments::print(out);
    }
    void usage();  // pretty printing of values/defaults/etc
    void getoptH();
};

class RemapArguments : public BaseArguments
{
 public:
    std::string chromosome;
    int insertSize;
    int occurrenceCutoff;
    int wordSize;

    int readLength;
    int readSkipOffset;

    std::string                 outputFilename;
    std::vector<std::string>    references;
    std::vector<std::string>    SAMHeaderOptions;
 RemapArguments() :
    insertSize(250),
        occurrenceCutoff(5000),
        wordSize(15),
        readLength(100),
        readSkipOffset(100)
        {
            optionString += "a:c:l:o:O:r:s:w:";
        }
    ~RemapArguments()
    {
        ;
    }
    void print(std::ostream &out)
    {
        out << "chromosome " << chromosome << " (defaults to all)" << std::endl;
        out << "insertSize " << insertSize << " (default 250)" << std::endl;
        out << "occurrenceCutoff " << occurrenceCutoff << " (default 5000)" << std::endl;
        out << "reference (no default):";
        for (std::vector<std::string>::iterator it=references.begin(); it<references.end(); it++)
        {
            out << " " << *it;
        }
        out << "readLength " << readLength << " (default 100)" << std::endl;
        out << "readSkipOffset " << readSkipOffset << " (default 100)" << std::endl;
        out << std::endl;

        BaseArguments::print(out);
    }
    void usage();  // pretty printing of values/defaults/etc
    void getoptH();
};

class CreateArguments : public BaseArguments
{
 public:
    bool createIndex;
    bool isColorSpace;
    int occurrenceCutoff;
    int wordSize;
 CreateArguments() :
    createIndex(false),
        isColorSpace(false),
        occurrenceCutoff(5000),
        wordSize(15)
        {
            optionString += "ciO:w:";
        }
    ~CreateArguments()
    {
        ;
    }
    void print(std::ostream &out)
    {
        out << "createIndex " << (createIndex ? "on" : "off") << " (default off)" << std::endl;
        out << "isColorSpace " << (isColorSpace ? "on" : "off") << " (default off)" << std::endl;
        out << "occurrenceCutoff " << occurrenceCutoff << " (default 5000)" << std::endl;
        out << "wordSize " << wordSize << " (default 15)" << std::endl;
        BaseArguments::print(out);
    }
    void usage();  // pretty printing of values/defaults/etc
    void getoptH();
};

class CheckArguments : public BaseArguments
{
 public:
    bool verbose;
 CheckArguments() : verbose(false)
    {
        optionString += "v";
    }
    ~CheckArguments()
    {
        ;
    }
    void print(std::ostream &out)
    {
        out << "verbose " << (verbose ? "on" : "off") << " (default off)" << std::endl;
        BaseArguments::print(out);
    }
    void usage();  // pretty printing of values/defaults/etc
    void getoptH();
};

class HeaderArguments : public BaseArguments
{
 public:
    bool editFlag;
    bool isColorSpace;
    std::string newHeaderFile;
 HeaderArguments() : editFlag(false), isColorSpace(false)
    {
        optionString += "ceh:";
    }
    ~HeaderArguments()
    {
        ;
    }
    void print(std::ostream &out)
    {
        out << "editFlag " << (editFlag ? "on" : "off") << " (default off)" << std::endl;
        out << "isColorSpace " << (isColorSpace ? "on" : "off") << " (default off)" << std::endl;
        out << "newHeaderFile " << newHeaderFile << std::endl;
        BaseArguments::print(out);
    }
    void usage();  // pretty printing of values/defaults/etc
    void getoptH();
};

//
// Unused at the moment.
//
class FastQCheckArguments : public BaseArguments
{
 public:
    bool isColorSpace;
    int minReadLength;
    int maxReportedErrors;
 FastQCheckArguments() : isColorSpace(false), minReadLength(10), maxReportedErrors(20)
    {
        optionString += "m:e:c";
    }
    ~FastQCheckArguments()
    {
        ;
    }
    void print(std::ostream &out)
    {
        out << "baseSpace? " << (!isColorSpace? "base space":"color space") << " (default base space)"<<std::endl;
        out << "minReadLength " << minReadLength << " (default 10)" << std::endl;
        out << "maxReportErrors " << maxReportedErrors << " (default 20)" << std::endl;
        BaseArguments::print(out);
    }
    void usage();  // pretty printing of values/defaults/etc
    void getoptH();
};

#endif
