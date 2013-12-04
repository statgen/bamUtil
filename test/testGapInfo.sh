#!/bin/bash

status=0;
../bin/bam gapInfo --in testFiles/testGapInfo.sam --out results/gapInfo.txt --noph 2> results/gapInfo.log
let "status |= $?"
diff results/gapInfo.txt expected/gapInfo.txt
let "status |= $?"
diff results/gapInfo.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testGapInfo.sam --refFile testFilesLibBam/chr1_partial.fa --out results/gapInfoRef.txt --noph 2> results/gapInfoRef.log
let "status |= $?"
diff results/gapInfoRef.txt expected/gapInfoRef.txt
let "status |= $?"
diff results/gapInfoRef.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testGapInfo.sam --detailed --out results/gapInfoDetailed.txt --noph 2> results/gapInfoDetailed.log
let "status |= $?"
diff results/gapInfoDetailed.txt expected/gapInfoDetailed.txt
let "status |= $?"
diff results/gapInfoDetailed.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testGapInfo.sam --detailed --checkFirst --checkStrand --out results/gapInfoDetailedCheck.txt --noph 2> results/gapInfoDetailedCheck.log
let "status |= $?"
diff results/gapInfoDetailedCheck.txt expected/gapInfoDetailedCheck.txt
let "status |= $?"
diff results/gapInfoDetailedCheck.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testGapInfo.sam --detailed --checkFirst --out results/gapInfoDetailedCheckFirst.txt --noph 2> results/gapInfoDetailedCheckFirst.log
let "status |= $?"
diff results/gapInfoDetailedCheckFirst.txt expected/gapInfoDetailedCheckFirst.txt
let "status |= $?"
diff results/gapInfoDetailedCheckFirst.log expected/empty.log
let "status |= $?"

../bin/bam gapInfo --in testFiles/testGapInfo.sam --detailed --checkStrand --out results/gapInfoDetailedCheckStrand.txt --noph 2> results/gapInfoDetailedCheckStrand.log
let "status |= $?"
diff results/gapInfoDetailedCheckStrand.txt expected/gapInfoDetailedCheckStrand.txt
let "status |= $?"
diff results/gapInfoDetailedCheckStrand.log expected/empty.log
let "status |= $?"


if [ $status != 0 ]
then
  echo failed testGapInfo.sh
  exit 1
fi

