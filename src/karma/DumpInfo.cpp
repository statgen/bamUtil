
#include "DumpInfo.h"
#include "FastQFile.h"
#include "GenomeSequence.h"
#include <string>

#define HAS_SUFFIX(filename,suffix) \
    (strlen(suffix) < strlen(filename) && \
     strcasecmp(filename + strlen(filename) - strlen(suffix), suffix)==0)

void DumpInfo(CheckArguments args, std::ostream &stream, const char *file)
{
    if (HAS_SUFFIX(file, ".bam"))
    {
        stream << "bam file!!\n";
    }
    else if (HAS_SUFFIX(file, ".sam"))
    {
    }
    else if (HAS_SUFFIX(file, ".fastq"))
    {
        BaseAsciiMap::SPACE_TYPE myBaseType = BaseAsciiMap::UNKNOWN;
        String filename = file;
        FastQFile validator;
        validator.validateFastQFile(filename, true, myBaseType);
    }
    else if (HAS_SUFFIX(file, ".umfa"))
    {
        GenomeSequence gs;
        if (gs.open(file))
        {
            stream << "Sequence file " << file << " could not be opened." << std::endl;
            stream << gs.getErrorString() << std::endl;
        }
        else
        {
            stream << gs.getHeader();
            gs.close();
        }
    }
    else if (HAS_SUFFIX(file, ".fa"))
    {
    }
    else if (HAS_SUFFIX(file, ".umwhl"))
    {
#if 0
        std::cout << std::endl << "Dumping Left Word Hash file header: " << std::endl;
        std::cout << whLeft.getHeader() << std::endl;
#endif
    }
    else if (HAS_SUFFIX(file, ".umwhr"))
    {
#if 0
        std::cout << std::endl << "Dumping Left Word Hash file header: " << std::endl;
        std::cout << whLeft.getHeader() << std::endl;
#endif
    }
    else if (HAS_SUFFIX(file, ".umwiwp"))
    {
#if 0
        std::cout << std::endl << "Dumping Word Index file header: " << std::endl;
        std::cout << wi << std::endl;
#endif

    }
    else if (HAS_SUFFIX(file, ".umwihi"))
    {
    }
    else
    {
    }
}


