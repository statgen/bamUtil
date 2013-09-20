#!/bin/bash

status=0;
../bin/bam dedup --in testFiles/testDedup.sam --out results/testDedup.sam 2> results/testDedup.txt
let "status |= $?"
diff results/testDedup.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup.sam expected/testDedup.sam
let "status |= $?"
diff results/testDedup.sam.log expected/testDedup.sam.log
let "status |= $?"

../bin/bam dedup --recab --in testFiles/testDedup.sam --out results/testDedupRecab.sam --refFile testFiles/ref_partial.fa 2> results/testDedupRecab.txt
let "status |= $?"
diff results/testDedupRecab.txt expected/testDedupRecab.txt
let "status |= $?"
diff results/testDedupRecab.sam expected/testDedupRecab.sam
let "status |= $?"
diff results/testDedupRecab.sam.log expected/testDedupRecab.sam.log
let "status |= $?"
sort results/testDedupRecab.sam.qemp | diff - expected/testDedupRecab.sam.qemp
let "status |= $?"

../bin/bam recab --in results/testDedup.sam --out results/testDedupRecab2Step.sam --refFile testFiles/ref_partial.fa 2> results/testDedupRecab2Step.txt
let "status |= $?"
diff results/testDedupRecab2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedupRecab2Step.sam expected/testDedupRecab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedupRecab2Step.sam.log expected/testDedupRecab2Step.sam.log
let "status |= $?"
sort results/testDedupRecab2Step.sam.qemp | diff - expected/testDedupRecab.sam.qemp
let "status |= $?"

../bin/bam dedup --in testFiles/testDedup1.sam --out results/testDedup1.sam 2> results/testDedup1.txt
let "status |= $?"
diff results/testDedup1.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup1.sam expected/testDedup1.sam
let "status |= $?"
diff results/testDedup1.sam.log expected/testDedup1.sam.log
let "status |= $?"

../bin/bam dedup --recab --in testFiles/testDedup1.sam --out results/testDedup1Recab.sam --refFile testFiles/ref_partial.fa 2> results/testDedup1Recab.txt
let "status |= $?"
diff results/testDedup1Recab.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup1Recab.sam expected/testDedup1Recab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup1Recab.sam.log expected/testDedup1Recab.sam.log
let "status |= $?"
sort results/testDedup1Recab.sam.qemp | diff - expected/testDedup1Recab.sam.qemp
let "status |= $?"

../bin/bam recab --in results/testDedup1.sam --out results/testDedup1Recab2Step.sam --refFile testFiles/ref_partial.fa 2> results/testDedup1Recab2Step.txt
let "status |= $?"
diff results/testDedup1Recab2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedup1Recab2Step.sam expected/testDedup1Recab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup1Recab2Step.sam.log expected/testDedup1Recab2Step.sam.log
let "status |= $?"
sort results/testDedup1Recab2Step.sam.qemp | diff - expected/testDedup1Recab.sam.qemp
let "status |= $?"

if [ $status != 0 ]
then
  echo failed testDedup.sh
  exit 1
fi

