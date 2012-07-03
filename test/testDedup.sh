#!/bin/bash

status=0;
../bin/bam dedup --in testFiles/testDedup.sam --out results/testDedup.sam 2> results/testDedup.txt
let "status |= $?"
diff results/testDedup.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup.sam expected/testDedup.sam
let "status |= $?"

../bin/bam recab --in testFiles/testRecab.sam --out results/testRecab.sam --refFile testFilesLibBam/chr1_partial.fa > results/testRecab.txt 2> results/testRecab.log
let "status |= $?"
diff results/testRecab.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecab.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecab.log expected/empty.log
let "status |= $?"



if [ $status != 0 ]
then
  echo failed testDedup.sh
  exit 1
fi

