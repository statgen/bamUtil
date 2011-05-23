# Diff identical sam/bam
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff1.bam --seq --baseQual --tags "OP:i;MD:Z" > results/diffSameSamBam.log 2> results/empty.log && diff results/diffSameSamBam.log expected/diffSameSamBam.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSam.log 2> results/empty.log && diff results/diffOrderSam.log expected/diffOrderSam.log && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh.log 2> results/empty.log && diff results/diffNoThresh.log expected/diffNoThresh.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --out results/diffOrderSam2.log 2> results/empty.log && diff results/diffOrderSam2.log expected/diffOrderSam2.log && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh2.log 2> results/empty.log && diff results/diffNoThresh2.log expected/diffNoThresh2.log && diff results/empty.log expected/empty.txt  \
&& \
################# Output as bam/sam.
# Diff identical sam/bam
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff1.bam --seq --baseQual --tags "OP:i;MD:Z" --out results/diffSameSamBam.bam 2> results/empty.log \
&& diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --onlyDiffs --out results/diffOrderSam.bam 2> results/empty.log \
&& diff results/diffOrderSam.bam expected/diffOrderSam.bam && diff results/diffOrderSam_only2_testDiff2.bam expected/diffOrderSam_only2_testDiff2.bam && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh.bam 2> results/empty.log \
&& diff results/diffNoThresh_only1_testDiff1.bam expected/diffNoThresh_only1_testDiff1.bam && diff results/diffNoThresh_only2_testDiff2.bam expected/diffNoThresh_only2_testDiff2.bam && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --out results/diffOrderSam2.sam 2> results/empty.log \
&& diff results/diffOrderSam2.sam expected/diffOrderSam2.sam && diff results/diffOrderSam2_only1_testDiff2.sam expected/diffOrderSam2_only1_testDiff2.sam && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos reverse the order.
../bin/bam diff --in1 testFiles/testDiff2.sam --in2 testFiles/testDiff1.sam --seq --baseQual --tags "OP:i;MD:Z" --posDiff 1 --out results/diffNoThresh2.sam 2> results/empty.log \
&& diff results/diffNoThresh2_only1_testDiff2.sam expected/diffNoThresh2_only1_testDiff2.sam && diff results/diffNoThresh2_only2_testDiff1.sam expected/diffNoThresh2_only2_testDiff1.sam && diff results/empty.log expected/empty.txt 

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
