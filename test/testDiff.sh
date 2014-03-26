# Diff identical sam/bam
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff1.bam --seq --baseQual --tags "OP:i;MD:Z" --noph > results/diffSameSamBam.log 2> results/empty.log && diff results/diffSameSamBam.log expected/diffSameSamBam.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSam.log --noph 2> results/empty.log && diff results/diffOrderSam.log expected/diffOrderSam.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records, but only 4 records in the pool.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 4 --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSamPool4.txt --noph 2> results/diffOrderSamPool4.log && diff results/diffOrderSamPool4.txt expected/diffOrderSamPool4.txt && diff results/diffOrderSamPool4.log expected/diffOrderSamPool4.log \
&& \
# Different order/pos on one of the records, but only 3 records in the pool.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 3 --seq --baseQual --tags "OP:i,MD:Z" --onlyDiffs --out results/diffOrderSamPool3.txt --noph 2> results/diffOrderSamPool3.log && diff results/diffOrderSamPool3.txt expected/diffOrderSamPool3.txt && diff results/diffOrderSamPool3.log expected/diffOrderSamPool3.log \
&& \
# Different order/pos on one of the records, but only 2 records in the pool.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 2 --seq --baseQual --tags "OP:i,MD:Z" --onlyDiffs --out results/diffOrderSamPool2.txt --noph 2> results/diffOrderSamPool2.log && diff results/diffOrderSamPool2.txt expected/diffOrderSamPool2.txt && diff results/diffOrderSamPool2.log expected/diffOrderSamPool2.log \
&& \
# Diff order/pos but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh.log --noph 2> results/empty.log && diff results/diffNoThresh.log expected/diffNoThresh.log && diff results/empty.log expected/empty.txt \
&& \
# Diff all but no threshold for pos.
../bin/bam diff --all --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i,MD:Z" --posDiff 1 --out results/diffNoThreshAll.log --noph 2> results/empty.log && diff results/diffNoThreshAll.log expected/diffNoThreshAll.log && diff results/empty.log expected/empty.txt \
&& \
# Diff all but no threshold for pos, diffs only.
../bin/bam diff --all --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --posDiff 1 --onlyDiffs --out results/diffNoThreshAllDiffOnly.log --noph 2> results/empty.log && diff results/diffNoThreshAllDiffOnly.log expected/diffNoThreshAllDiffOnly.log && diff results/empty.log expected/empty.txt \
&& \
# Diff flag only but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --flag --noCigar --noPos --posDiff 1 --out results/diffFlag.log --noph 2> results/empty.log && diff results/diffFlag.log expected/diffFlag.log && diff results/empty.log expected/empty.txt \
&& \
# Diff mapqual, mate, and isize only but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --mapQual --mate --isize --noCigar --noPos --posDiff 1 --out results/diffMateMapQISize.log --noph 2> results/empty.log && diff results/diffMateMapQISize.log expected/diffMateMapQISize.log && diff results/empty.log expected/empty.txt \
&& \
# Diff only all tags but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --everyTag --noCigar --noPos --posDiff 1 --out results/diffTagsOnly.log --noph 2> results/empty.log && diff results/diffTagsOnly.log expected/diffTagsOnly.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --out results/diffOrderSam2.log --noph 2> results/empty.log && diff results/diffOrderSam2.log expected/diffOrderSam2.log && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh2.log --noph 2> results/empty.log && diff results/diffNoThresh2.log expected/diffNoThresh2.log && diff results/empty.log expected/empty.txt  \
&& \
################# Output as bam/sam.
# Diff identical sam/bam
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff1.bam --seq --baseQual --tags "OP:i;MD:Z" --out results/diffSameSamBam.bam --noph 2> results/empty.log \
&& diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSam.bam --noph 2> results/empty.log \
&& diff results/diffOrderSam.bam expected/diffOrderSam.bam && diff results/diffOrderSam_only2_testDiff2.bam expected/diffOrderSam_only2_testDiff2.bam && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records but with poolsize of 4.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 4 --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSamPool4.sam --noph 2> results/diffOrderSamPool4.log \
&& diff results/diffOrderSamPool4.sam expected/diffOrderSamPool4.sam && diff results/diffOrderSamPool4_only2_testDiff2.sam expected/diffOrderSamPool4_only2_testDiff2.sam && diff results/diffOrderSamPool4.log expected/diffOrderSamPool4.log \
&& \
# Different order/pos on one of the records but with poolsize of 3.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 3 --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSamPool3.sam --noph 2> results/diffOrderSamPool3.log \
&& diff results/diffOrderSamPool3.sam expected/diffOrderSamPool3.sam && diff results/diffOrderSamPool3_only2_testDiff2.sam expected/diffOrderSamPool3_only2_testDiff2.sam && diff results/diffOrderSamPool3.log expected/diffOrderSamPool3.log \
&& \
# Different order/pos on one of the records but with poolsize of 2.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 2 --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSamPool2.sam --noph 2> results/diffOrderSamPool2.log \
&& diff results/diffOrderSamPool2_only1_testDiff1.sam expected/diffOrderSamPool2_only1_testDiff1.sam && diff results/diffOrderSamPool2_only2_testDiff2.sam expected/diffOrderSamPool2_only2_testDiff2.sam && diff results/diffOrderSamPool2.log expected/diffOrderSamPool2.log \
&& \
# Different order/pos on one of the records but with poolsize of 1.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --recPoolSize 1 --seq --baseQual --tags "OP:i,MD:Z" --onlyDiffs --out results/diffOrderSamPool1.sam --noph 2> results/diffOrderSamPool1.log \
&& diff results/diffOrderSamPool1.sam expected/diffOrderSamPool1.sam && diff results/diffOrderSamPool1_only1_testDiff1.sam expected/diffOrderSamPool1_only1_testDiff1.sam && diff results/diffOrderSamPool1_only2_testDiff2.sam expected/diffOrderSamPool1_only2_testDiff2.sam && diff results/diffOrderSamPool1.log expected/diffOrderSamPool1.log \
&& \
# Diff order/pos but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh.bam --noph 2> results/empty.log \
&& diff results/diffNoThresh_only1_testDiff1.bam expected/diffNoThresh_only1_testDiff1.bam && diff results/diffNoThresh_only2_testDiff2.bam expected/diffNoThresh_only2_testDiff2.bam && diff results/empty.log expected/empty.txt \
&& \
# Diff all but no threshold for pos.
../bin/bam diff --all --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThreshAll.sam --noph 2> results/empty.log && diff results/diffNoThreshAll.sam expected/diffNoThreshAll.sam && diff results/diffNoThreshAll_only1_testDiff1.sam expected/diffNoThreshAll_only1_testDiff1.sam && diff results/diffNoThreshAll_only2_testDiff2.sam expected/diffNoThreshAll_only2_testDiff2.sam && diff results/empty.log expected/empty.txt \
&& \
# Diff all but no threshold for pos, diffs only.
../bin/bam diff --all --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --posDiff 1 --onlyDiffs --out results/diffNoThreshAllDiffOnly.log --noph 2> results/empty.log && diff results/diffNoThreshAllDiffOnly.log expected/diffNoThreshAllDiffOnly.log && diff results/empty.log expected/empty.txt \
&& \
# Diff flag only but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --flag --noCigar --noPos --posDiff 1 --out results/diffFlag.log --noph 2> results/empty.log && diff results/diffFlag.log expected/diffFlag.log && diff results/empty.log expected/empty.txt \
&& \
# Diff mapqual, mate, and isize only but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --mapQual --mate --isize --noCigar --noPos --posDiff 1 --out results/diffMateMapQISize.log --noph 2> results/empty.log && diff results/diffMateMapQISize.log expected/diffMateMapQISize.log && diff results/empty.log expected/empty.txt \
&& \
# Diff only all tags but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --everyTag --noCigar --noPos --posDiff 1 --out results/diffTagsOnly.log --noph 2> results/empty.log && diff results/diffTagsOnly.log expected/diffTagsOnly.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --out results/diffOrderSam2.sam --noph 2> results/empty.log \
&& diff results/diffOrderSam2.sam expected/diffOrderSam2.sam && diff results/diffOrderSam2_only1_testDiff2.sam expected/diffOrderSam2_only1_testDiff2.sam && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh2.sam --noph 2> results/empty.log \
&& diff results/diffNoThresh2_only1_testDiff2.sam expected/diffNoThresh2_only1_testDiff2.sam && diff results/diffNoThresh2_only2_testDiff1.sam expected/diffNoThresh2_only2_testDiff1.sam && diff results/empty.log expected/empty.txt \
&& \
################ Test for ASHG poster 2011
# Diff semi-real example diff output.
../bin/bam diff --in1 testFiles/testDiff3.sam --in2 testFiles/testDiff4.sam --seq --baseQual --tags "NM:i;OQ:Z;OP:i;OC:Z" --out results/diffAshg.log --noph 2> results/empty.log \
&& diff results/diffAshg.log expected/diffAshg.log && diff results/empty.log expected/empty.txt \
#&& \
# Diff semi-real example sam output.
../bin/bam diff --in1 testFiles/testDiff3.sam --in2 testFiles/testDiff4.sam --seq --baseQual --tags "NM:i;OQ:Z;OP:i;OC:Z" --out results/diffAshg.sam --noph 2> results/empty.log \
&& diff results/diffAshg.sam expected/diffAshg.sam && diff results/diffAshg_only1_testDiff3.sam expected/diffAshg_only1_testDiff3.sam && diff results/diffAshg_only2_testDiff4.sam expected/diffAshg_only2_testDiff4.sam && diff results/empty.log expected/empty.txt \

if [ $? -ne 0 ]
then
  exit 1
fi

# Check for the existence of files that should not be there
if [ -e results/diffSameSamBam.bam ]
then
  exit 1
fi
if [ -e results/diffSameSamBam_only1_testDiff1.bam ]
then
  exit 1
fi
if [ -e results/diffSameSamBam_only2_testDiff1.bam ]
then
  exit 1
fi
if [ -e results/diffOrderSam_only1_testDiff1.bam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool4_only1_testDiff1.sam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool3_only1_testDiff1.sam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool2.sam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool1.sam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool1_only1_testDiff1.sam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool1_only2_testDiff2.sam ]
then
  exit 1
fi
if [ -e results/diffOrderSamPool1.sam ]
then
  exit 1
fi
if [ -e results/diffNoThresh.bam ]
then
  exit 1
fi
if [ -e results/diffOrderSam2_only2_testDiff1.sam ]
then
  exit 1
fi
if [ -e results/diffNoThresh2.sam ]
then
  exit 1
fi
