#!/bin/bash

#####
# The qemp files must be softed prior to diffing because the code that generates
# them loops through either a map or an unordered map.  The order of the
# unordered map is undefined and varies by implementation of it.


status=0;
###############
# Recalibration
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecab.sam --refFile testFilesLibBam/chr1_partial.fa --fitModel > results/testRecab.txt 2> results/testRecab.log
let "status |= $?"
diff results/testRecab.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecab.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecab.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecab.sam.qemp) <(sort expected/testRecab.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecab.sam.log expected/testRecab.sam.log
let "status |= $?"

../bin/bam recab --noph --fast --in testFiles/testRecab.sam --out results/testRecabFast.sam --refFile testFilesLibBam/chr1_partial.fa --fitModel > results/testRecabFast.txt 2> results/testRecabFast.log
let "status |= $?"
diff results/testRecabFast.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecabFast.txt expected/empty.txt
let "status |= $?"
diff results/testRecabFast.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabFast.sam.qemp) <(sort expected/testRecabFast.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabFast.sam.log expected/testRecabFast.sam.log
let "status |= $?"

# Store the original quality
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabStoreQ.sam --refFile testFilesLibBam/chr1_partial.fa --storeQualTag OQ --fitModel > results/testRecabStoreQ.txt 2> results/testRecabStoreQ.log
let "status |= $?"
diff results/testRecabStoreQ.sam expected/testRecabStoreQ.sam
let "status |= $?"
diff results/testRecabStoreQ.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecabStoreQ.log expected/testRecabStoreQ.log
let "status |= $?"
diff <(sort results/testRecabStoreQ.sam.qemp) <(sort expected/testRecab.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabStoreQ.sam.log expected/testRecabStoreQ.sam.log
let "status |= $?"

# Use OQ
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabUseOQ.sam --refFile testFilesLibBam/chr1_partial.fa --qualField OQ --fitModel > results/testRecabUseOQ.txt 2> results/testRecabUseOQ.log
let "status |= $?"
diff results/testRecabUseOQ.sam expected/testRecabUseOQ.sam
let "status |= $?"
diff results/testRecabUseOQ.txt expected/testRecabUseOQ.txt
let "status |= $?"
diff results/testRecabUseOQ.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabUseOQ.sam.qemp) <(sort expected/testRecabUseOQ.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabUseOQ.sam.log expected/testRecabUseOQ.sam.log
let "status |= $?"

# Use OQ & store
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabUseStoreOQ.sam --refFile testFilesLibBam/chr1_partial.fa --qualField OQ --storeQualTag OQ > results/testRecabUseStoreOQ.txt --fitModel 2> results/testRecabUseStoreOQ.log
let "status |= $?"
diff results/testRecabUseStoreOQ.sam expected/testRecabUseStoreOQ.sam
let "status |= $?"
diff results/testRecabUseStoreOQ.txt expected/testRecabUseOQ.txt
let "status |= $?"
diff results/testRecabUseStoreOQ.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabUseStoreOQ.sam.qemp) <(sort expected/testRecabUseOQ.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabUseStoreOQ.sam.log expected/testRecabUseStoreOQ.sam.log
let "status |= $?"

# Bad quality tag
../bin/bam recab --noph --in testFiles/testRecabBad.sam --out results/testRecabUseBadOQ.sam --refFile testFilesLibBam/chr1_partial.fa --qualField OQ --storeQualTag OQ --fitModel > results/testRecabUseBadOQ.txt 2> results/testRecabUseBadOQ.log
let "status |= $?"
diff results/testRecabUseBadOQ.sam expected/testRecabUseBadOQ.sam
let "status |= $?"
diff results/testRecabUseBadOQ.txt expected/testRecabUseBadOQ.txt
let "status |= $?"
diff results/testRecabUseBadOQ.log expected/testRecabUseBadOQ.log
let "status |= $?"
diff <(sort results/testRecabUseBadOQ.sam.qemp) <(sort expected/testRecabUseBadOQ.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabUseBadOQ.sam.log expected/testRecabUseBadOQ.sam.log
let "status |= $?"
#TODO

# Quality tag does not exist
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabNoTag.sam --refFile testFilesLibBam/chr1_partial.fa --qualField ZQ --fitModel > results/testRecabNoTag.txt 2> results/testRecabNoTag.log
let "status |= $?"
diff results/testRecabNoTag.sam expected/testRecab.sam
let "status |= $?"
diff results/testRecabNoTag.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecabNoTag.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabNoTag.sam.qemp) <(sort expected/testRecab.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabNoTag.sam.log expected/testRecabNoTag.sam.log
let "status |= $?"


###############
# Test with DBSNP
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabDBSNP.sam --refFile testFilesLibBam/chr1_partial.fa --dbsnp testFiles/dbsnp1.txt > results/testRecabDBSNP.txt 2> results/testRecabDBSNP.log
let "status |= $?"
diff results/testRecabDBSNP.sam expected/testRecabDBSNP.sam
let "status |= $?"
diff results/testRecabDBSNP.txt expected/empty.txt
let "status |= $?"
diff results/testRecabDBSNP.log expected/testRecabDBSNP.log
let "status |= $?"
diff <(sort results/testRecabDBSNP.sam.qemp) <(sort expected/testRecabDBSNP.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabDBSNP.sam.log expected/testRecabDBSNP.sam.log
let "status |= $?"

###############
# Test with DBSNP.gz
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabDBSNPgz.sam --refFile testFilesLibBam/chr1_partial.fa --dbsnp testFiles/dbsnp1.txt.gz > results/testRecabDBSNPgz.txt 2> results/testRecabDBSNPgz.log
let "status |= $?"
diff results/testRecabDBSNPgz.sam expected/testRecabDBSNP.sam
let "status |= $?"
diff results/testRecabDBSNPgz.txt expected/empty.txt
let "status |= $?"
diff results/testRecabDBSNPgz.log expected/testRecabDBSNPgz.log
let "status |= $?"
diff <(sort results/testRecabDBSNPgz.sam.qemp) <(sort expected/testRecabDBSNP.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabDBSNPgz.sam.log expected/testRecabDBSNPgz.sam.log
let "status |= $?"

###############
# Test with DBSNP, keeping even if previous is dbsnp
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabDBSNPkeepPrev.sam --refFile testFilesLibBam/chr1_partial.fa --dbsnp testFiles/dbsnp1.txt --keepPrevDbsnp > results/testRecabDBSNPkeepPrev.txt 2> results/testRecabDBSNPkeepPrev.log
let "status |= $?"
diff results/testRecabDBSNPkeepPrev.sam expected/testRecabDBSNPkeepPrev.sam
let "status |= $?"
diff results/testRecabDBSNPkeepPrev.txt expected/empty.txt
let "status |= $?"
diff results/testRecabDBSNPkeepPrev.log expected/testRecabDBSNPkeepPrev.log
let "status |= $?"
diff <(sort results/testRecabDBSNPkeepPrev.sam.qemp) <(sort expected/testRecabDBSNPkeepPrev.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabDBSNPkeepPrev.sam.log expected/testRecabDBSNPkeepPrev.sam.log
let "status |= $?"

###############
# Recalibration with indels
../bin/bam recab --noph --in testFiles/testRecab2.sam --out results/testRecab2.sam --refFile testFilesLibBam/chr1_partial.fa > results/testRecab2.txt 2> results/testRecab2.log
let "status |= $?"
diff results/testRecab2.sam expected/testRecab2.sam
let "status |= $?"
diff results/testRecab2.txt expected/empty.txt
let "status |= $?"
diff results/testRecab2.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecab2.sam.qemp) <(sort expected/testRecab2.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecab2.sam.log expected/testRecab2.sam.log
let "status |= $?"

###############
# Recalibration with indels keeping even if previous is not match/mismatch
../bin/bam recab --noph --in testFiles/testRecab2.sam --out results/testRecab2Keep.sam --keepPrevNonAdjacent --refFile testFilesLibBam/chr1_partial.fa > results/testRecab2Keep.txt 2> results/testRecab2Keep.log
let "status |= $?"
diff results/testRecab2Keep.sam expected/testRecab2Keep.sam
let "status |= $?"
diff results/testRecab2Keep.txt expected/empty.txt
let "status |= $?"
diff results/testRecab2Keep.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecab2Keep.sam.qemp) <(sort expected/testRecab2Keep.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecab2Keep.sam.log expected/testRecab2Keep.sam.log
let "status |= $?"


###############
# Recalibration max quality
../bin/bam recab --noph --in testFiles/testRecab51.sam --out results/testRecabMaxQual50.sam --refFile testFilesLibBam/chr1_partial.fa > results/testRecabMaxQual50.txt 2> results/testRecabMaxQual50.log
let "status |= $?"
diff results/testRecabMaxQual50.sam expected/testRecabMaxQual50.sam
let "status |= $?"
diff results/testRecabMaxQual50.txt expected/empty.txt
let "status |= $?"
diff results/testRecabMaxQual50.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabMaxQual50.sam.qemp) <(sort expected/testRecabMaxQual51.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabMaxQual50.sam.log expected/testRecabMaxQual50.sam.log
let "status |= $?"

###############
# Recalibration min quality
../bin/bam recab --noph --in testFiles/testRecab51.sam --out results/testRecabMinQual28.sam --minBaseQual 28 --refFile testFilesLibBam/chr1_partial.fa > results/testRecabMinQual28.txt 2> results/testRecabMinQual28.log
let "status |= $?"
diff results/testRecabMinQual28.sam expected/testRecabMinQual28.sam
let "status |= $?"
diff results/testRecabMinQual28.txt expected/empty.txt
let "status |= $?"
diff results/testRecabMinQual28.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabMinQual28.sam.qemp) <(sort expected/testRecabMinQual28.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabMinQual28.sam.log expected/testRecabMinQual28.sam.log
let "status |= $?"

###############
# Recalibration max quality
../bin/bam recab --noph --in testFiles/testRecab51.sam --out results/testRecabMaxQual51.sam --maxBaseQual 51 --refFile testFilesLibBam/chr1_partial.fa > results/testRecabMaxQual51.txt 2> results/testRecabMaxQual51.log
let "status |= $?"
diff results/testRecabMaxQual51.sam expected/testRecabMaxQual51.sam
let "status |= $?"
diff results/testRecabMaxQual51.txt expected/empty.txt
let "status |= $?"
diff results/testRecabMaxQual51.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabMaxQual51.sam.qemp) <(sort expected/testRecabMaxQual51.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabMaxQual51.sam.log expected/testRecabMaxQual51.sam.log
let "status |= $?"

###############
# Recalibration max quality
../bin/bam recab --noph --in testFiles/testRecab51.sam --out results/testRecabMaxQual49.sam --maxBaseQual 49 --refFile testFilesLibBam/chr1_partial.fa > results/testRecabMaxQual49.txt 2> results/testRecabMaxQual49.log
let "status |= $?"
diff results/testRecabMaxQual49.sam expected/testRecabMaxQual49.sam
let "status |= $?"
diff results/testRecabMaxQual49.txt expected/empty.txt
let "status |= $?"
diff results/testRecabMaxQual49.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabMaxQual49.sam.qemp) <(sort expected/testRecabMaxQual51.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabMaxQual49.sam.log expected/testRecabMaxQual49.sam.log
let "status |= $?"

###############
# Recalibration with binning
../bin/bam recab --noph --in testFiles/testRecab.sam --out results/testRecabBin.sam --refFile testFilesLibBam/chr1_partial.fa --fitModel --binMid --binQualS 2,3,10,20,25,30,35,40,50 > results/testRecabBin.txt 2> results/testRecabBin.log
let "status |= $?"
diff results/testRecabBin.sam expected/testRecabBin.sam
let "status |= $?"
diff results/testRecabBin.txt expected/testRecab.txt
let "status |= $?"
diff results/testRecabBin.log expected/empty.log
let "status |= $?"
diff <(sort results/testRecabBin.sam.qemp) <(sort expected/testRecab.sam.qemp)
let "status |= $?"
diff -I "Start: .*" -I "End: .*" results/testRecabBin.sam.log expected/testRecabBin.sam.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed testRecab.sh
  exit 1
fi

