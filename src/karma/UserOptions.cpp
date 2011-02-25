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

#include "UserOptions.h"
#include <stdint.h>
#include <unistd.h>

//
// Handle arguments that are common to all internal command modes.
//
// Derived classes need to call the base class getopt in order
// to ensure that all arguments are handled.  Base class issues
// errors for unrecognized arguments.
//
void BaseArguments::getopt()
{
    optind = 1;
    while ((opt = ::getopt(argc, (char **)(argv), (const char *) optionString.c_str())) != -1)
    {
        //
        // Go down the class hierarchy, checking options.
        // First one to get the option wins, and returns.
        //
        getoptH();  // go down the class hierarchy, checking options
    }
}

void BaseArguments::getoptH()
{
    switch (opt)
    {
        case 'd':
            debug = true;
            return;
        case 's':
            seed = atoi(optarg);
            return;
        default:
            usage();
            exit(1);
            // XXX NOTREACHED
    }
}

void BaseArguments::usage()
{
    std::cerr  << "   -d           -> debug" << std::endl;
    std::cerr  << "   -s [int]     -> set random number seed [12345]" << std::endl;

    std::cerr  << std::endl;
    std::cerr  << "Defaults:" << std::endl;
    std::cerr  << std::endl;
    this->print(std::cerr);
    exit(1);
}

////////////////////////////////////////////////////////////
//
// "Map" argument methods
//
void MapArguments::usage()
{
    std::cerr << "Usage:" << std::endl;
    std::cerr  << "   " << argv[0] << " [options...] file1.fastq.gz [file2.fastq.gz]" << std::endl;
    std::cerr  << "   -a [int]     -> maximum insert size" << std::endl;
    std::cerr  << "   -B [int]     -> max number of bases in millions" << std::endl;
    std::cerr  << "   -c           -> alignment in color space (ABI SOLiD data)" << std::endl;
    std::cerr  << "   -E           -> show reference bases (default off)" << std::endl;
    std::cerr  << "   -H           -> set SAM header line values (e.g. -H RG:SM:NA12345)" << std::endl;
    std::cerr  << "   -o           -> required output sam/bam file" << std::endl;
    std::cerr  << "   -O [int]     -> occurrence cutoff (default 5000)" << std::endl;
    std::cerr  << "   -q [int]     -> soft clip read bases with low quality" << std::endl;
    std::cerr  << "   -Q           -> quiet mode (no output except errors)" << std::endl;
    std::cerr  << "   -r [name]    -> required genome reference" << std::endl;
    std::cerr  << "   -R [int]     -> max number of reads" << std::endl;
    std::cerr  << "   -w [int]     -> index word size (default 15)" << std::endl;
    std::cerr  << std::endl;
    std::cerr  << "-O and -w modify which reference index gets opened." << std::endl;
    std::cerr  << std::endl;
    std::cerr  << "Simplest examples:" << std::endl;
    std::cerr  << std::endl;
    std::cerr  << "karma map -r NCBI36.fa -o result.sam mate1.fastq.gz mate2.fastq.gz" << std::endl;
    std::cerr  << "karma map -r NCBI36.fa -o single.sam single_end.fastq.gz" << std::endl;
    std::cerr  << "karma map -r NCBI36.fa -o - mate1.fastq.gz mate2.fastq.gz" << std::endl;
    BaseArguments::usage();
}

void MapArguments::getoptH()
{
    switch (opt)
    {
        case 'a':
            insertSize = atoi(optarg);
            return;
        case 'B':
            maxBases = atoi(optarg) * 1000 * 1000;
            return;
        case 'c':
            mapInColorSpace = true;
            return;
        case 'E':
            showReferenceBases = true;
            return;
        case 'H':
            SAMHeaderOptions.push_back(optarg);
            return;
        case 'o':
            outputFilename = optarg;
            return;
        case 'O':
            occurrenceCutoff = atoi(optarg);
            return;
        case 'q':
            qualityTrim = atoi(optarg);
            return;
        case 'Q':
            quietMode = true;
            return;
        case 'r':
            references.push_back(optarg);
            return;
        case 'R':
            maxReads = atoi(optarg);
            return;
        case 'w':
            wordSize = atoi(optarg);
            return;
            break;
    }
    BaseArguments::getoptH();
}

//
// "Remap" argument methods
//
void RemapArguments::usage()
{
    std::cerr << "Usage:" << std::endl;
    std::cerr  << "   " << argv[0] << " [options...]" << std::endl;
    std::cerr  << "   -a [int]     -> maximum insert size" << std::endl;
    std::cerr  << "   -c           -> optional chromosome to process (default all)" << std::endl;
    std::cerr  << "   -o           -> optional output file" << std::endl;
    std::cerr  << "   -O [int]     -> occurrence cutoff (default 5000)" << std::endl;
    std::cerr  << "   -r [name]    -> required genome reference" << std::endl;
    std::cerr  << "   -s [int]     -> read skip offset (default 100)" << std::endl;
    std::cerr  << std::endl;
    BaseArguments::usage();
}

void RemapArguments::getoptH()
{
    switch (opt)
    {
        case 'a':
            insertSize = atoi(optarg);
            return;
        case 'c':
            chromosome = optarg;
            return;
        case 'o':
            outputFilename = optarg;
            return;
        case 'O':
            occurrenceCutoff = atoi(optarg);
            return;
        case 'r':
            references.push_back(optarg);
            return;
        case 's':
            readSkipOffset = atoi(optarg);
            return;
        case 'w':
            wordSize = atoi(optarg);
            return;
            break;
    }
    BaseArguments::getoptH();
}

//
// "Create" argument methods
//
void CreateArguments::usage()
{
    std::cerr << "Usage:" << std::endl;
    std::cerr  << "   " << argv[0] << " [options...] reference.fa" << std::endl;
    std::cerr  << "   -c           -> isColorSpace" << std::endl;
    std::cerr  << "   -i           -> create index used for map command" << std::endl;
    std::cerr  << "   -w           -> # of bases in word size" << std::endl;
    BaseArguments::usage();
}

void CreateArguments::getoptH()
{
    switch (opt)
    {
        case 'c':
            isColorSpace = true;
            return;
        case 'i':
            createIndex = true;
            return;
        case 'w':
            wordSize = atoi(optarg);
            return;
            break;
    }
    BaseArguments::getoptH();
}

//
// "Check" argument methods
//
void CheckArguments::usage()
{
    std::cerr << "Usage:" << std::endl;
    std::cerr  << "   " << argv[0] << " [options...] file.{sam|bam|fa|fastq|umfa} [...]" << std::endl;
    std::cerr  << "   -v           -> verbose" << std::endl;
    BaseArguments::usage();
}

void CheckArguments::getoptH()
{
    switch (opt)
    {
        case 'v':
            verbose = true;
            return;
            break;
    }
    BaseArguments::getoptH();
}

//
// "Header" argument methods
//
void HeaderArguments::usage()
{
    std::cerr << "Usage:" << std::endl;
    std::cerr  << "   " << argv[0] << " [options...] file.{fa|umfa}" << std::endl;
    std::cerr  << "   -c           -> operate on the color space reference" << std::endl;
    std::cerr  << "   -e           -> edit header" << std::endl;
    std::cerr  << "   -h <filename>  -> load new header from this file" << std::endl;
    BaseArguments::usage();
}

void HeaderArguments::getoptH()
{
    switch (opt)
    {
        case 'c':
            isColorSpace = true;
            return;
        case 'e':
            editFlag = true;
            return;
        case 'h':
            newHeaderFile = optarg;
            return;
            break;
    }
    BaseArguments::getoptH();
}

//
// "FastQCheck" argument methods
//
void FastQCheckArguments::usage()
{
    std::cerr << "Usage:" << std::endl;
    std::cerr  << "   " << argv[0] << " [options...] file " << std::endl;
    std::cerr  << "   -c -> set to color space" << std::endl;
    std::cerr  << "   -m -> mininum read length" << std::endl;
    std::cerr  << "   -e -> maximum output errors" << std::endl;
    BaseArguments::usage();
}

void FastQCheckArguments::getoptH()
{
    switch (opt)
    {
        case 'b':
            isColorSpace = false;
            return;
        case 'm':
            minReadLength = atoi(optarg);
            return;
        case 'e':
            maxReportedErrors = atoi(optarg);
            return;
    }
    BaseArguments::getoptH();
}
