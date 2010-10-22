#ifndef __FASTQ_STATUS_H__
#define __FASTQ_STATUS_H__

#include <string>

class FastQStatus
{
 public:

   // Return value enum for the FastQFile class methods.
   // Enum to indicate the error codes.
   //   FASTQ_SUCCESS indicates method finished successfully.
   //   FASTQ_INVALID means that the sequence was invalid.
   //   FASTQ_ORDER_ERROR means the methods are called out of order, like trying to read a file before opening it.
   //   FASTQ_OPEN_ERROR means the file could not be opened.
   //   FASTQ_CLOSE_ERROR means the file could not be closed.
   //   FASTQ_READ_ERROR means that a problem occurred on a read.
   //   FASTQ_NO_SEQUENCE_ERROR means there were no errors, but no sequences read.
   enum Status {FASTQ_SUCCESS = 0, FASTQ_INVALID, FASTQ_ORDER_ERROR, FASTQ_OPEN_ERROR, FASTQ_CLOSE_ERROR, FASTQ_READ_ERROR, FASTQ_NO_SEQUENCE_ERROR};

   // Get the enum string for the status.
   static const char* getStatusString(Status status);

private:
   static const char* enumString[];
};


#endif
