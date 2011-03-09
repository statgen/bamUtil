#ifndef _FASTQREADER_H_
#define _FASTQREADER_H_

#include "InputFile.h"
#include "Error.h"
#include "StringBasics.h"
struct Fastq {
    String tag;
    String read;
    String tag2;
    String qual;

    bool isValidFastqRecord() {
        if (this->tag.Length() <= 1 || this->tag[0] != '@')
            error("Tag line must start with [@]"); //. check line %lu\n", this->lineNo);
        if (this->read.Length() == 0)
            error("unexpected empty short read DNA strand");// line %lu\n", this->lineNo);
        if (this->qual.Length() == 0)
            error("unexpected empty DNA quality");// line %lu\n", this->lineNo);
        return true;
    };

};

class FastqReader{
 private:
    IFILE fHandle;
 public:
    FastqReader(const char* fileName) {
        this->fHandle = ifopen(fileName, "rb");
        if (!this->fHandle)
            error("Cannot open file: %s", fileName);
    };
    void Close() {
        ifclose(this->fHandle);
    };
    bool Eof() {
        return ifeof(this->fHandle);
    }

 public:
    // Read from f atmost nReads fastq reads into buffer
    // @param IFILE f: an openned IFILE handle
    // @param buffer Fastq* buffer: non-null memory buffer
    // @param int nReads: at most this number of reads will be read into memory
    // @return int: how many Fastq reads are actully read, 0: readed end; minus values: error!
    int ReadFastqFile(Fastq* buffer, int nReads);
};

#endif /* _FASTQREADER_H_ */
