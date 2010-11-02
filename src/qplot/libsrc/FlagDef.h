#ifndef __FLAGDEF_H__
#define __FLAGDEF_H__

//BAM_FPAIRED
//
//the read is paired in sequencing, no matter whether it is mapped in a pair
//
#define BAM_FPAIRED 1 
//
//BAM_FPROPER_PAIR
//
//the read is mapped in a proper pair
//
#define BAM_FPROPER_PAIR 2 
//
//BAM_FUNMAP
//
//the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
//
#define BAM_FUNMAP 4 
//
//BAM_FMUNMAP
//
//the mate is unmapped
//
#define BAM_FMUNMAP 8 
//
//BAM_FREVERSE
//
//the read is mapped to the reverse strand
//
#define BAM_FREVERSE 16 
//
//BAM_FMREVERSE
//
//the mate is mapped to the reverse strand
//
#define BAM_FMREVERSE 32 
//
//BAM_FREAD1
//
//this is read1
//
#define BAM_FREAD1 64 
//
//BAM_FREAD2
//
//this is read2
//
#define BAM_FREAD2 128 
//
//BAM_FSECONDARY
//
//not primary alignment
//
#define BAM_FSECONDARY 256 
//
//BAM_FQCFAIL
//
//QC failure
//
#define BAM_FQCFAIL 512 
//
//BAM_FDUP
//
//optical or PCR duplicate
//
#define BAM_FDUP 1024 

#endif
