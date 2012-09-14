#!/bin/bash

status=0;
###############
# Recalibration
../bin/bam recab --in testFiles/testRecab.sam --out results/testRecab.sam --refFile testFilesLibBam/chr1_partial.fa > results/testRecab.txt 2> results/testRecab.log
let "status |= $?"
diff results/testRecab.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecab.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecab.log expected/empty.log
let "status |= $?"
#diff results/testRecab.sam.qemp expected/testRecab.sam.qemp
#let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecab.sam.log expected/testRecab.sam.log
let "status |= $?"

../bin/bam recab --fast --in testFiles/testRecab.sam --out results/testRecabFast.sam --refFile testFilesLibBam/chr1_partial.fa > results/testRecabFast.txt 2> results/testRecabFast.log
let "status |= $?"
diff results/testRecabFast.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecabFast.txt expected/empty.txt
let "status |= $?"
diff results/testRecabFast.log expected/empty.log
let "status |= $?"
#diff results/testRecab.sam.qemp expected/testRecab.sam.qemp
#let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabFast.sam.log expected/testRecabFast.sam.log
let "status |= $?"

# Store the original quality
../bin/bam recab --in testFiles/testRecab.sam --out results/testRecabStoreQ.sam --refFile testFilesLibBam/chr1_partial.fa --storeQualTag OQ > results/testRecabStoreQ.txt 2> results/testRecabStoreQ.log
let "status |= $?"
diff results/testRecabStoreQ.sam expected/testRecabStoreQ.sam
let "status |= $?"
diff results/testRecabStoreQ.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecabStoreQ.log expected/empty.log
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabStoreQ.sam.log expected/testRecabStoreQ.sam.log
let "status |= $?"

# Use OQ
../bin/bam recab --in testFiles/testRecab.sam --out results/testRecabUseOQ.sam --refFile testFilesLibBam/chr1_partial.fa --qualField OQ > results/testRecabUseOQ.txt 2> results/testRecabUseOQ.log
let "status |= $?"
diff results/testRecabUseOQ.sam expected/testRecabUseOQ.sam
let "status |= $?"
diff results/testRecabUseOQ.txt expected/testRecabUseOQ.txt
let "status |= $?"
diff results/testRecabUseOQ.log expected/empty.log
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabUseOQ.sam.log expected/testRecabUseOQ.sam.log
let "status |= $?"

# Use OQ & store
../bin/bam recab --in testFiles/testRecab.sam --out results/testRecabUseStoreOQ.sam --refFile testFilesLibBam/chr1_partial.fa --qualField OQ --storeQualTag OQ > results/testRecabUseStoreOQ.txt 2> results/testRecabUseStoreOQ.log
let "status |= $?"
diff results/testRecabUseStoreOQ.sam expected/testRecabUseStoreOQ.sam
let "status |= $?"
diff results/testRecabUseStoreOQ.txt expected/testRecabUseOQ.txt
let "status |= $?"
diff results/testRecabUseStoreOQ.log expected/empty.log
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabUseStoreOQ.sam.log expected/testRecabUseStoreOQ.sam.log
let "status |= $?"

# Bad quality tag
../bin/bam recab --in testFiles/testRecabBad.sam --out results/testRecabUseBadOQ.sam --refFile testFilesLibBam/chr1_partial.fa --qualField OQ --storeQualTag OQ > results/testRecabUseBadOQ.txt 2> results/testRecabUseBadOQ.log
let "status |= $?"
diff results/testRecabUseBadOQ.sam expected/testRecabUseBadOQ.sam
let "status |= $?"
diff results/testRecabUseBadOQ.txt expected/testRecabUseBadOQ.txt
let "status |= $?"
diff results/testRecabUseBadOQ.log expected/empty.log
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabUseBadOQ.sam.log expected/testRecabUseBadOQ.sam.log
let "status |= $?"
#TODO

# Quality tag does not exist
../bin/bam recab --in testFiles/testRecab.sam --out results/testRecabNoTag.sam --refFile testFilesLibBam/chr1_partial.fa --qualField ZQ > results/testRecabNoTag.txt 2> results/testRecabNoTag.log
let "status |= $?"
diff results/testRecabNoTag.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecabNoTag.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecabNoTag.log expected/empty.log
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabNoTag.sam.log expected/testRecabNoTag.sam.log
let "status |= $?"


###############
# Test with DBSNP
../bin/bam recab --in testFiles/testRecab.sam --out results/testRecabDBSNP.sam --refFile testFilesLibBam/chr1_partial.fa --dbsnp testFiles/dbsnp1.txt --skipFit > results/testRecabDBSNP.txt 2> results/testRecabDBSNP.log
let "status |= $?"
diff results/testRecabDBSNP.sam expected/testRecabDBSNP.sam
let "status |= $?"
diff results/testRecabDBSNP.txt expected/empty.txt
let "status |= $?"
diff results/testRecabDBSNP.log expected/testRecabDBSNP.log
let "status |= $?"
#diff results/testRecabDBSNP.sam.qemp expected/testRecab.sam.qemp
#let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabDBSNP.sam.log expected/testRecabDBSNP.sam.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed testRecab.sh
  exit 1
fi

