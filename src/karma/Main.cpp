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

#include "DumpInfo.h"
#include "Error.h"
#include "GenomeSequence.h"
#include "MapperPE.h"
#include "MapperPEBaseSpace.h"
#include "ReadsProcessor.h"
#include "Test.h"
#include "Random.h"
#include "Util.h"
#include "WordIndex.h"
#include "SamHeader.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <unistd.h>

#include "Main.h"

//
// somewhat messy routine to get the filenames using the basename
// obtained from the appropriate reference, and the arguments wordSize
// and occurrenceCutoff.
//
static std::string getReferenceNameWithArgs(
    int wordSize,
    int occurrenceCutoff,
    GenomeSequence *baseSpaceReference,
    GenomeSequence *colorSpaceReference)
{
    std::ostringstream buf;
    // ugly, but we want to open the color space index if it was provided,
    // otherwise the base space one:
    GenomeSequence &referenceTmp = colorSpaceReference ? *colorSpaceReference : *baseSpaceReference;

    buf << referenceTmp.getBaseFilename()
    << "-" << (referenceTmp.isColorSpace() ? "cs" : "bs") << "."
    << wordSize << "." << occurrenceCutoff;
    return buf.str();
}


void mainCheck(const char *program, int argc, const char **argv)
{
    CheckArguments args;
    args.setArgs(argc, argv);
    args.getopt();

    if (argc < 2)
    {
        args.usage();   // no return
    }

    //
    // given various files of various types, do various checks
    // *.umfa - print version, date, who made it, etc
    // *.fa - # chromosomes, lengths, etc.
    // *.umwhwp - size of word pointers table, etc
    // *.umwhhi - table loading?  other info
    // *.sam, *.bam - basic info?
    // *.fastq - similar?
    //

    for (int i=optind; i<argc; i++)
    {
        DumpInfo(args, std::cout, argv[i]);
    }
}

/**
 * Obtain a usable header
 */
static int parseHeader(SamHeader& header, MapArguments& args, std::string& commandLine) {
    header.clear();

    // get other fields for the SAM output header:
    header.set("HD", "VN", "1.0");
    header.set("HD", "SO", "unsorted");
    header.set("PG", "ID", "karma");
    header.set("PG", "VN", VERSION);
    header.set("PG", "CL", commandLine.c_str());

    // bool isReadGroupSet = false; // if RG tag is not set, we will set it to a default value.

    //
    // now handle SAM header overrides from the command line:
    //
    for (std::vector<std::string>::iterator it=args.SAMHeaderOptions.begin(); it<args.SAMHeaderOptions.end(); it++)
    {
        // format is for example RG:SM:NA12345
        // note that since DT can have colons in it, we only pay attention to the first two colons
        int firstColon = (*it).find_first_of(':');
        if (firstColon!=2)
        {
            std::cerr << "Bad SAM header option: " << *it << std::endl;
            args.usage();   // no return
        }

        int secondColon = (*it).find_first_of(':', firstColon+1);
        if (secondColon!=5)
        {
            std::cerr << "Bad SAM header option: " << *it << std::endl;
            args.usage();   // no return
        }
        std::string headerLineName = (*it).substr(0, firstColon);   // should always be length 2
        std::string headerLineTag = (*it).substr(firstColon+1, secondColon - firstColon - 1);   // should always be length 2
        std::string headerLineValue = (*it).substr(secondColon+1);   // arbitrary length

        if (headerLineValue.size()==0)
        {
            std::cerr << "Bad SAM header option (value size is zero): " << *it << std::endl;
            args.usage();   // no return
        }

        // // ugly way to track this... not sure of best solution.
        // if (headerLineName=="RG" && headerLineTag=="ID")
        //     isReadGroupSet = true;

        header.set(headerLineName, headerLineTag, headerLineValue);
    }


    // validate SAMHeader:
    // check if starred tag exists (to conform to SAM specification)
    if (!header.conformSpecification())
        return 1;


    // check if RG tag is present, since...
    // get the date for the SAM header:
    if (!header.containTag("RG")) {
        time_t t = time(NULL);
        struct tm *timep = gmtime(&t);
        char timeBuffer[256];
        strftime(timeBuffer,sizeof(timeBuffer), "%Y-%m-%dT%H:%MZ", timep);

        //
        // RG is special - if the line exists, both SM and ID must also
        // exist, even if set to unknown.  I decided I always wanted
        // DT set, so I also have to set SM and ID here.
        //
        //
        const char* defaultReadGroupID = "1";
        header.set("RG", "DT", timeBuffer);
        header.set("RG", "SM", "Unknown");
        header.set("RG", "ID", defaultReadGroupID);
    }

    // engine.mapperOptions.showReferenceBases = true;
    return 0;
}


void mainMap(const char *program, int argc, const char **argv)
{
    MapArguments args;
    args.setArgs(argc, argv);
    args.getopt();

    if (args.outputFilename == "-") args.quietMode = true;

    // we expect one or two sequence files
    if (optind >= argc || optind < argc - 2)
    {
        args.usage();   // no return
    }

    if (args.outputFilename=="")
    {
        std::cerr << "Must specify an output filename using -o." << std::endl;
        args.usage();
    }

    // only allow user to specify one reference name, since we now have unified way to open base/color space ref. genome.
    if (args.references.size() != 1)
    {
        std::cerr << "Unspecified base space reference" << std::endl;
        args.usage();
        exit(1);
    }

    ReadsProcessor engine;
    // get the commandline string for the SAM header and stdout:
    std::string commandLine = program;
    for (int i=0; i<argc; i++)
    {
        commandLine += " ";
        commandLine+=argv[i];
    }
    SamHeader header; 
    if ( parseHeader(header, args, commandLine) ) { 
        std::cerr << "Cannot set SAM header." << std::endl;
    }
    engine.setHeader(header);

    // setting mapper related parameters
    engine.parseMapArguments(args);

    // open reference, wordindex, wordhash
    if (engine.openReference(args.references[0], args.wordSize, args.occurrenceCutoff, args.quietMode, args.debug)) {
        std::cerr << "Open reference failed! " << std::endl;
        exit(1);
    }

    std::string sequenceFilename1 = argv[optind];
    std::string sequenceFilename2;

    if (optind == argc - 2)
    {
        sequenceFilename2 = argv[optind+1];
    }

    if (sequenceFilename2!="")
    {
        // paired end mapping
        engine.MapPEReadsFromFiles(
            sequenceFilename1,
            sequenceFilename2,
            args.outputFilename);
    }
    else
    {
        // single end mapping
        engine.MapSEReadsFromFile(
            sequenceFilename1,
            args.outputFilename);
    }
    engine.closeReference();
}

void mainCreate(const char *program, int argc, const char **argv)
{
    CreateArguments args;
    args.setArgs(argc, argv);

    args.getopt();

    if (optind >= argc)
    {
        args.usage();   // no return
    }

    std::string fileName = argv[::optind];

    GenomeSequence reference;

    reference.setProgressStream(std::cout);

    // create the reference UMFA file
    if (reference.setReferenceName(fileName))
    {
        std::cerr << "failed to set the genome reference file " << fileName << "." << std::endl;
        exit(1);
    }

    reference.setApplication("KARMA");

    //
    // once the reference is created with isColorSpace, the
    // value percolates through the wordIndex and then ReadsProcessor,
    // so it is sufficient to set it here on creation of the reference only.
    //
    reference.setColorSpace(args.isColorSpace);
    if (reference.create())
    {
        // error
        std::cerr << "failed to create the reference file " << fileName << "." << std::endl;
        exit(1);
    }

#if 0
    //
    // XXX A LOT OF THINGS NEED TO BE FILLED IN HERE:
    //
    // the header now exists, and is read/write, so we set
    // these up here (putting them before ::create() above,
    // these initializations will fail).
    //
    reference.setAssemblyID(args.assemblyID);
    reference.setSpecies(args.species);
    reference.setURI(args.uri);
#endif

    if (args.createIndex)
    {
        //
        // WordIndex creates two files, WordHash creates
        // a single file (one for each direction).
        //
        WordIndex wi;
        WordHash whLeft;
        WordHash whRight;

        //
        // Pass in a few paramaters for WordIndex creation
        //
        wi.setWordSize(args.wordSize);
        wi.setOccurrenceCutoff(args.occurrenceCutoff);
        wi.setVerbose(false);
        wi.setDebug(args.debug);

        //
        // Create the word index and word pointers table.
        // This is fairly slow, and also creates a map
        // in order that whLeft/whRight can created, so
        // is quite large.
        //
        std::cout << "creating reference short word index files." << std::endl;

        std::string    baseName = getReferenceNameWithArgs(args.wordSize, args.occurrenceCutoff, &reference, NULL);
        wi.setFilenames(baseName.c_str());

        wi.create(reference, whLeft, whRight, args.wordSize);

        //
        // Create the left hash.  The left hash is a
        // hash from a (2 * wordSize)-mer word to the possible
        // locations where it can be mapped.  The hash
        // key is created from the "wordSize" bases to the
        // left of the index word which has more than
        // occurrenceCutoff positions in the genome plus
        // the "wordsize" bases in the index word itself.
        //
        // The right hash is similar to the above left hash.
        //
        // This part of the process could be made smaller by
        // having WordIndex close the word index and word
        // position tables, I think.
        //
        std::string leftHashName, rightHashName;

        leftHashName = baseName + ".umwhl";
        rightHashName = baseName + ".umwhr";

        std::cout << "Creating left hash:\n" << std::endl;
        whLeft.create(leftHashName.c_str());

        std::cout << "Creating right hash:\n" << std::endl;
        whRight.create(rightHashName.c_str());
    }
}

void mainHeader(const char *program, int argc, const char **argv)
{
    HeaderArguments args;
    args.setArgs(argc, argv);

    args.getopt();

    // allow 1 additional argument (the reference)
    if (optind >= argc || optind < argc - 1)
    {
        args.usage();   // no return
    }

    std::string referenceFilename = argv[::optind];

    GenomeSequence reference;

    reference.setReferenceName(referenceFilename);
    reference.open(args.isColorSpace, (args.editFlag ? O_RDWR: O_RDONLY));
    if (args.editFlag)
    {
        // first write the TSV data to a temporary file
        char tmpFilename[TMP_MAX];
        strcpy(tmpFilename, "/tmp/karma_header.XXXXXX");
        int result = mkstemp(tmpFilename);
        if (result<0)
        {
            perror("mkstemp");
            exit(1);
        }
        ::close(result);
        std::ofstream editFile;
        editFile.open(tmpFilename);
        reference.dumpHeaderTSV(editFile);
        editFile.close();
        // now fire off $(EDITOR) or vi
        std::string editor = getenv("EDITOR");
        std::string editCommand;
        if (editor=="") editor = "vi";
        editCommand = editor + " " + tmpFilename;
        system(editCommand.c_str());
        // now ask if we really want to update the header
        std::cout << "Are you sure you want to update the header? [Ny] ";
        std::string yesNo;
        std::cin >> yesNo;
        if (yesNo=="yes" || yesNo=="Yes")
        {
            // read the header file, update the header
        }
        else
        {
            // save the edit file?
            // std::cout << "Your edits are stored in " << tmpFileName;
        }
    }
    else
    {
        reference.dumpHeaderTSV(std::cout);
    }
}

void mainRemap(const char *program, int argc, const char **argv)
{
    RemapArguments args;
    args.setArgs(argc, argv);
    args.getopt();

    // optind points to remainder

    GenomeSequence reference;
    reference.setDebugFlag(args.debug);

    WordIndex wi;

    wi.setVerbose(false);
    wi.setDebug(args.debug);

    if (args.references.size()!=1)
    {
        args.usage();
    }

    reference.setReferenceName(args.references[0].c_str());

    std::cerr << "open and prefetch reference: " << std::flush;

    // we always need to open the base space reference:
    if (reference.open(false))
    {
        std::cerr << "failed.\n";
        exit(1);
    }

    std::string    baseName = getReferenceNameWithArgs(args.wordSize, args.occurrenceCutoff, &reference, NULL);
    wi.setFilenames(baseName.c_str());

    std::cerr << "open and prefetch index and position tables: " << std::flush;

    if (wi.open(reference))
    {
        std::cerr << "failed.\n";
        exit(1);
    }

    std::cerr << std::endl;

    std::cerr << "open and prefetch long word hashes: " << std::flush;

    WordHash whLeft;
    WordHash whRight;

    std::string leftHashName, rightHashName;

    leftHashName = baseName + ".umwhl";
    rightHashName = baseName + ".umwhl";

    // open the word hashes
    if (whLeft.open(leftHashName.c_str()))
    {
        std::cerr << "failed to open left word hash " << leftHashName << "." << std::endl;
        exit(1);
    }

    if (whRight.open(rightHashName.c_str()))
    {
        std::cerr << "failed to open right word hash " << rightHashName  << "." << std::endl;
    }

    std::cerr << std::endl;

    std::cerr << "Now generating remapped reads from genome." << std::endl << std::flush;

    Test test;
    test.setGenomeSequence(&reference);
    test.setWordIndex(&wi);
    test.setWordHashLeft(&whLeft);
    test.setWordHashRight(&whRight);

    test.testRemapReference(args.outputFilename, args.chromosome, args.readLength, args.readSkipOffset);

}

void mainTest(const char *program, int argc, const char **argv)
{
    MapArguments args;
    args.setArgs(argc, argv);
    args.getopt();

    // optind points to remainder

    GenomeSequence reference;
    reference.setDebugFlag(args.debug);

    WordIndex wi;

    wi.setVerbose(false);
    wi.setDebug(false);

    reference.setReferenceName("/data/local/ref/karma.ref/human_b36_male.fa");

    // we always need to open the base space reference:
    if (reference.open(false))
    {
        std::cerr << "failed.\n";
        exit(1);
    }

    std::string    baseName = getReferenceNameWithArgs(args.wordSize, args.occurrenceCutoff, &reference, NULL);
    wi.setFilenames(baseName.c_str());

    if (wi.open(reference))
    {
        std::cerr << "failed.\n";
        exit(1);
    }

    std::cout << "open and prefetch long word hashes: " << std::flush;
    WordHash whLeft;
    WordHash whRight;

    std::string leftHashName, rightHashName;

    leftHashName = baseName + ".umwhl";
    rightHashName = baseName + ".umwhl";

    // open the word hashes
    if (whLeft.open(leftHashName.c_str()))
    {
        std::cerr << "failed to open left word hash " << leftHashName << "." << std::endl;
        exit(1);
    }
    if (whRight.open(rightHashName.c_str()))
    {
        std::cerr << "failed to open right word hash " << rightHashName  << "." << std::endl;
    }

    Test test;
    test.setGenomeSequence(&reference);
    test.setWordIndex(&wi);
    test.setWordHashLeft(&whLeft);
    test.setWordHashRight(&whRight);

    test.test(); //test base space, SE & PE reads.

}

void commandList(std::ostream &out, int argc, const char **argv)
{
    out << "Usage:  " << std::endl;
    out << "   " << argv[0] << " create [options...]" << std::endl;
    out << "   " << argv[0] << " header [options...] <reference genome>.umfa" << std::endl;
    out << "   " << argv[0] << " map [options...] file1.fastq.gz [file2.fastq.gz]" << std::endl;
    out << "Diagnostics:" << std::endl;
#if 0
    out << "   " << argv[0] << " check [options...]" << std::endl;
#endif
    out << "   " << argv[0] << " test [options...]" << std::endl;
    out << "   " << argv[0] << " remap [options...]" << std::endl;
    out << std::endl;
    out << "To get help for a Karma sub-command, just type e.g. 'karma map'" << std::endl;
}

int main(int argc, const char ** argv)
{
    BaseArguments args;
    args.setArgs(argc, argv);

    if (argc<2)
    {
        std::cerr << "Welcome to K-tuple Alignment with Rapid Matching Algorithm version " << VERSION << "." << std::endl;
        std::cout << "Written by Paul Anderson, Xiaowei Zhan, and Yun Li at The University of Michigan." << std::endl;
        std::cerr << "Built on " << DATE << " on host " << NODE << " by " << USER << "." << std::endl;
        std::cerr << std::endl;
        std::cerr << "Karma documentation, support and downloads are at: http://genome.sph.umich.edu/wiki/Karma" << std::endl;
        std::cerr << std::endl;
        commandList(std::cerr, argc, argv);
        args.usage();   // no return
    }

    if (strcmp(argv[1], "map")==0) mainMap(argv[0], argc - 1, argv + 1);
    else if (strcmp(argv[1], "remap")==0) mainRemap(argv[0], argc - 1, argv + 1);
#if 0
    else if (strcmp(argv[1], "check")==0) mainCheck(argv[0], argc - 1, argv + 1);
#endif
    else if (strcmp(argv[1], "create")==0) mainCreate(argv[0], argc - 1, argv + 1);
    else if (strcmp(argv[1], "header")==0) mainHeader(argv[0], argc - 1, argv + 1);
    else if (strcmp(argv[1], "test")==0) mainTest(argv[0], argc - 1, argv + 1);
    else
    {
        std::cerr << "Karma sub-command '" << argv[1] << "' not recognized." << std::endl;
        std::cerr << std::endl;
        commandList(std::cerr, argc, argv);
        args.usage();   // no return
    }
}
