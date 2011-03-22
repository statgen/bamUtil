#!/bin/sh

REFERENCE_GENOME=/home/sequences/build36
REFERENCE_DESTINATION=/home/1000G/data
REFERENCE_BASENAME=chromosomesTest
# fasta format file (ascii):
REFERENCE_DESTINATION_FA=$REFERENCE_DESTINATION/$REFERENCE_BASENAME.fa
# UM version of fasta (binary for memory mapping):
REFERENCE_DESTINATION_UMFA=$REFERENCE_DESTINATION/$REFERENCE_BASENAME.umfa

#
# Setting up the reference is done exactly once.
#
# First, a concatenated set of fasta files is created
# in ascending order.
#
# Second, karma is run to create the memory mapped binary
# version of the same file, plus word indices that are used
# during later mapping runs.
#
# Lastly, karma is run as many times as necessary to map the
# input read files.
#
# setup_reference should be safe to run multiple times - won't overwrite existing files:
#
setup_reference () {
	if [ ! -f $REFERENCE_DESTINATION_FA ] ; 
	then
		cat $REFERENCE_GENOME/Homo_sapiens.NCBI36.49.dna.chromosome.{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X}.fa.gz >$REFERENCE_DESTINATION_FA
	fi
	if [ ! -f $REFERENCE_DESTINATION_UMFA ] ;
	then
		./karma --create --reference $REFERENCE_DESTINATION_FA --readCutoff 5000
	fi
}


setup_reference

#
# test single end mapping of test data file.
#
./karma --reference $REFERENCE_DESTINATION_FA --maxReads 10000 ../testdata/test1.fastq.gz 

#
# test paired mapping of test data files
#
./karma --pairedReads --reference $REFERENCE_DESTINATION_FA --maxReads 10000 ../testdata/test_paired_1.fastq.gz  ../testdata/test_paired_2.fastq.gz

