ERROR=false

# Sam To Sam
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqNoOpt.sam && diff results/seqNoOpt.sam expected/seqOrig.sam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqOrig.sam --seqOrig && diff results/seqOrig.sam expected/seqOrig.sam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqBases.sam --seqBases && diff results/seqBases.sam expected/seqBases.sam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqEquals.sam --seqEquals && diff results/seqEquals.sam expected/seqEquals.sam \
&& \
# Bam to Sam
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqNoOptBam.sam && diff results/seqNoOptBam.sam expected/seqOrig.sam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqOrigBam.sam --seqOrig && diff results/seqOrigBam.sam expected/seqOrig.sam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqBasesBam.sam --seqBases && diff results/seqBasesBam.sam expected/seqBases.sam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqEqualsBam.sam --seqEquals && diff results/seqEqualsBam.sam expected/seqEquals.sam \
&& \
# Sam To Bam
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqNoOpt.bam && diff results/seqNoOpt.bam expected/seqOrig.bam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqOrig.bam --seqOrig && diff results/seqOrig.bam expected/seqOrig.bam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqBases.bam --seqBases && diff results/seqBases.bam expected/seqBases.bam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.sam --out results/seqEquals.bam --seqEquals && diff results/seqEquals.bam expected/seqEquals.bam \
&& \
# Bam To Bam
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqNoOptBam.bam && diff results/seqNoOptBam.bam expected/seqOrig.bam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqOrigBam.bam --seqOrig && diff results/seqOrigBam.bam expected/seqOrig.bam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqBasesBam.bam --seqBases && diff results/seqBasesBam.bam expected/seqBases.bam \
&& \
../../bin/bam convert --refFile testFiles/chr1_partial.fa --in testFiles/testFilter.bam --out results/seqEqualsBam.bam --seqEquals && diff results/seqEqualsBam.bam expected/seqEquals.bam

if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi

