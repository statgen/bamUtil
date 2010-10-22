#include "FastQStatus.h"

const char* FastQStatus::enumString[] = {"FASTQ_SUCCESS", "FASTQ_INVALID", "FASTQ_ORDER_ERROR", "FASTQ_OPEN_ERROR", "FASTQ_CLOSE_ERROR", "FASTQ_READ_ERROR", "FASTQ_NO_SEQUENCE_ERROR"};


const char* FastQStatus::getStatusString(Status status)
{
   return(enumString[status]);
}
