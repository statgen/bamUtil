#include "FastqReader.h"
int FastqReader::ReadFastqFile(Fastq* buffer, int nReads)
{
    int retCode = 0;
    int processedRead = 0;
    int lineNumber = 0;

    while (nReads) {
        if (ifeof(this->fHandle)) {
            return processedRead;
        }

        // read tag line in fastq file
        retCode = buffer[processedRead].tag.ReadLine(this->fHandle);
        if (retCode < 0 || ifeof(this->fHandle) ){
            error("Wrong tag in fastq file ");
        }
        lineNumber ++ ;

        // read 'read' line in fastq file
        retCode = buffer[processedRead].read.ReadLine(this->fHandle);
        if (retCode < 0 || ifeof(this->fHandle) ){
            error("Wrong read in fastq file ");
        }
        lineNumber ++ ;

        // read the 3rd line in fastq line
        retCode = buffer[processedRead].tag2.ReadLine(this->fHandle);
        if (retCode < 0 || ifeof(this->fHandle) ){
            error("Wrong tag2 in fastq file ");
        }
        lineNumber ++ ;

        // read the quality line in fastq line
        retCode = buffer[processedRead].qual.ReadLine(this->fHandle);
        if (retCode < 0 ){
            error("Wrong quality in fastq file ");
        }
        lineNumber ++ ;

        if (!(buffer[processedRead].isValidFastqRecord()))
            error("Invalid Fastq record between line %lu and line %lu", lineNumber-3, lineNumber);

        processedRead++;
        if (processedRead == nReads) break;
    }
    return processedRead;
}
