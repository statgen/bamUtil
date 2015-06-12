#!/bin/bash

status=0;
../bin/bam dedup --in testFiles/testDedup.sam --out results/testDedup.sam --noph 2> results/testDedup.txt
let "status |= $?"
diff results/testDedup.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup.sam expected/testDedup.sam
let "status |= $?"
diff results/testDedup.sam.log expected/testDedup.sam.log
let "status |= $?"

../bin/bam dedup --recab --in testFiles/testDedup.sam --out results/testDedupRecab.sam --refFile testFiles/ref_partial.fa --noph 2> results/testDedupRecab.txt
let "status |= $?"
diff results/testDedupRecab.txt expected/testDedupRecab.txt
let "status |= $?"
diff results/testDedupRecab.sam expected/testDedupRecab.sam
let "status |= $?"
diff results/testDedupRecab.sam.log expected/testDedupRecab.sam.log
let "status |= $?"
sort results/testDedupRecab.sam.qemp | diff - expected/testDedupRecab.sam.qemp
let "status |= $?"

../bin/bam recab --in results/testDedup.sam --out results/testDedupRecab2Step.sam --refFile testFiles/ref_partial.fa --noph 2> results/testDedupRecab2Step.txt
let "status |= $?"
diff results/testDedupRecab2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedupRecab2Step.sam expected/testDedupRecab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedupRecab2Step.sam.log expected/testDedupRecab2Step.sam.log
let "status |= $?"
sort results/testDedupRecab2Step.sam.qemp | diff - expected/testDedupRecab.sam.qemp
let "status |= $?"

../bin/bam dedup --in testFiles/testDedup1.sam --out results/testDedup1.sam --noph 2> results/testDedup1.txt
let "status |= $?"
diff results/testDedup1.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup1.sam expected/testDedup1.sam
let "status |= $?"
diff results/testDedup1.sam.log expected/testDedup1.sam.log
let "status |= $?"

../bin/bam dedup --recab --in testFiles/testDedup1.sam --out results/testDedup1Recab.sam --refFile testFiles/ref_partial.fa --noph 2> results/testDedup1Recab.txt
let "status |= $?"
diff results/testDedup1Recab.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup1Recab.sam expected/testDedup1Recab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup1Recab.sam.log expected/testDedup1Recab.sam.log
let "status |= $?"
sort results/testDedup1Recab.sam.qemp | diff - expected/testDedup1Recab.sam.qemp
let "status |= $?"

../bin/bam recab --in results/testDedup1.sam --out results/testDedup1Recab2Step.sam --refFile testFiles/ref_partial.fa --noph 2> results/testDedup1Recab2Step.txt
let "status |= $?"
diff results/testDedup1Recab2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedup1Recab2Step.sam expected/testDedup1Recab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup1Recab2Step.sam.log expected/testDedup1Recab2Step.sam.log
let "status |= $?"
sort results/testDedup1Recab2Step.sam.qemp | diff - expected/testDedup1Recab.sam.qemp
let "status |= $?"

../bin/bam dedup --in testFiles/testDedup2.sam --out results/testDedup2Force.sam --force --noph 2> results/testDedup2Force.txt
let "status |= $?"
diff results/testDedup2Force.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup2Force.sam expected/testDedup1.sam
let "status |= $?"
diff results/testDedup2Force.sam.log expected/testDedup2Force.sam.log
let "status |= $?"

../bin/bam dedup --in testFiles/testDedup2.sam --out results/testDedup2Exclude.sam --excludeFlags 0x704 --noph 2> results/testDedup2Exclude.txt
let "status |= $?"
diff results/testDedup2Exclude.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup2Exclude.sam expected/testDedup2Exclude.sam
let "status |= $?"
diff results/testDedup2Exclude.sam.log expected/testDedup2Exclude.sam.log
let "status |= $?"

../bin/bam recab --in results/testDedup2Force.sam --out results/testDedup2ForceRecab2Step.sam --refFile testFiles/ref_partial.fa --noph 2> results/testDedup2ForceRecab2Step.txt
let "status |= $?"
diff results/testDedup2ForceRecab2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedup2ForceRecab2Step.sam expected/testDedup1Recab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup2ForceRecab2Step.sam.log expected/testDedup2ForceRecab2Step.sam.log
let "status |= $?"
sort results/testDedup2ForceRecab2Step.sam.qemp | diff - expected/testDedup1Recab.sam.qemp
let "status |= $?"

../bin/bam dedup --recab --in testFiles/testDedup2.sam --force --out results/testDedup2ForceRecab.sam --refFile testFiles/ref_partial.fa --noph 2> results/testDedup2ForceRecab.txt
let "status |= $?"
diff results/testDedup2ForceRecab.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup2ForceRecab.sam expected/testDedup1Recab.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup2ForceRecab.sam.log expected/testDedup2ForceRecab.sam.log
let "status |= $?"
sort results/testDedup2ForceRecab.sam.qemp | diff - expected/testDedup1Recab.sam.qemp
let "status |= $?"

../bin/bam recab --in results/testDedup2Force.sam --out results/testDedup2ForceRecabD2Step.sam --refFile testFiles/ref_partial.fa --buildExclude 0x204 --noph 2> results/testDedup2ForceRecabD2Step.txt
let "status |= $?"
diff results/testDedup2ForceRecabD2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedup2ForceRecabD2Step.sam expected/testDedup2ForceRecabD.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup2ForceRecabD2Step.sam.log expected/testDedup2ForceRecabD2Step.sam.log
let "status |= $?"
sort results/testDedup2ForceRecabD2Step.sam.qemp | diff - expected/testDedup2ForceRecabD.sam.qemp
let "status |= $?"

../bin/bam dedup --recab --force --in testFiles/testDedup2.sam --out results/testDedup2ForceRecabD.sam --refFile testFiles/ref_partial.fa --buildExclude 0x204 --noph 2> results/testDedup2ForceRecabD.txt
let "status |= $?"
diff results/testDedup2ForceRecabD.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup2ForceRecabD.sam expected/testDedup2ForceRecabD.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup2ForceRecabD.sam.log expected/testDedup2ForceRecabD.sam.log
let "status |= $?"
sort results/testDedup2ForceRecabD.sam.qemp | diff - expected/testDedup2ForceRecabD.sam.qemp
let "status |= $?"

../bin/bam dedup --recab --force --in testFiles/testDedup2.sam --out results/testDedup2ForceRecabDA.sam --refFile testFiles/ref_partial.fa --excludeFlags 0x304 --buildExclude 0x204 --applyExclude 0x304 --noph 2> results/testDedup2ForceRecabDA.txt
let "status |= $?"
diff results/testDedup2ForceRecabDA.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup2ForceRecabDA.sam expected/testDedup2ForceRecabDA.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup2ForceRecabDA.sam.log expected/testDedup2ForceRecabDA.sam.log
let "status |= $?"
sort results/testDedup2ForceRecabDA.sam.qemp | diff - expected/testDedup2ForceRecabD.sam.qemp
let "status |= $?"

../bin/bam recab --in results/testDedup2Force.sam --out results/testDedup2ForceRecabDA2Step.sam --refFile testFiles/ref_partial.fa --buildExclude 0x204 --applyExclude 0x304 --noph 2> results/testDedup2ForceRecabDA2Step.txt
let "status |= $?"
diff results/testDedup2ForceRecabDA2Step.txt expected/empty.txt
let "status |= $?"
diff results/testDedup2ForceRecabDA2Step.sam expected/testDedup2ForceRecabDA.sam
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testDedup2ForceRecabDA2Step.sam.log expected/testDedup2ForceRecabDA2Step.sam.log
let "status |= $?"
sort results/testDedup2ForceRecabDA2Step.sam.qemp | diff - expected/testDedup2ForceRecabD.sam.qemp
let "status |= $?"

#Bin
../bin/bam dedup --recab --in testFiles/testDedup.sam --out results/testDedupRecabBin.sam --refFile testFiles/ref_partial.fa --noph --binMid --binQualS 2,3,10,20,25,30,35,40,50 2> results/testDedupRecabBin.txt
let "status |= $?"
diff results/testDedupRecabBin.txt expected/testDedupRecab.txt
let "status |= $?"
diff results/testDedupRecabBin.sam expected/testDedupRecabBin.sam
let "status |= $?"
diff results/testDedupRecabBin.sam.log expected/testDedupRecabBin.sam.log
let "status |= $?"
sort results/testDedupRecabBin.sam.qemp | diff - expected/testDedupRecab.sam.qemp
let "status |= $?"



if [ $status != 0 ]
then
  echo failed testDedup.sh
  exit 1
fi

