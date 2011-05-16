# Diff identical sam/bam
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff1.bam --baseQual > results/diffSameSamBam.log 2> results/empty.log && diff results/diffSameSamBam.log expected/diffSameSamBam.log && diff results/empty.log expected/empty.txt \
&& \
# Different order/pos on one of the records.
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --baseQual --out results/diffOrderSam.log 2> results/empty.log && diff results/diffOrderSam.log expected/diffOrderSam.log && diff results/empty.log expected/empty.txt \
&& \
# Diff order/pos but no threshold for pos
../bin/bam diff --in1 testFiles/testDiff1.sam --in2 testFiles/testDiff2.sam --baseQual --posDiff 1 --out results/diffNoThresh.log 2> results/empty.log && diff results/diffNoThresh.log expected/diffNoThresh.log && diff results/empty.log expected/empty.txt \

if [ $? -ne 0 ]
then
  exit 1
fi

