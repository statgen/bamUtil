#!/bin/bash

status=0;
../bin/bam gapInfo --in testFiles/testClipOverlapCoord.sam --out results/gapInfo.txt 2> results/gapInfo.log
let "status |= $?"
diff results/gapInfo.txt expected/gapInfo.txt
let "status |= $?"
diff results/gapInfo.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testClipOverlapCoord.sam --checkFirst --checkStrand --out results/gapInfoCheck.txt 2> results/gapInfoCheck.log
let "status |= $?"
diff results/gapInfoCheck.txt expected/gapInfoCheck.txt
let "status |= $?"
diff results/gapInfoCheck.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testClipOverlapCoord.sam --checkFirst --out results/gapInfoCheckFirst.txt 2> results/gapInfoCheckFirst.log
let "status |= $?"
diff results/gapInfoCheckFirst.txt expected/gapInfoCheckFirst.txt
let "status |= $?"
diff results/gapInfoCheckFirst.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testClipOverlapCoord.sam --checkStrand --out results/gapInfoCheckStrand.txt 2> results/gapInfoCheckStrand.log
let "status |= $?"
diff results/gapInfoCheckStrand.txt expected/gapInfoCheckStrand.txt
let "status |= $?"
diff results/gapInfoCheckStrand.log expected/empty.log
let "status |= $?"


if [ $status != 0 ]
then
  echo failed testGapInfo.sh
  exit 1
fi

