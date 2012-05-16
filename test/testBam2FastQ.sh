#!/bin/bash

status=0;
# Test clipping files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadName 2> results/testBam2FastQReadName.log
let "status |= $?"
diff results/testBam2FastQReadName.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQReadName_1.fastq expected/testBam2FastQReadName_1.fastq
let "status |= $?"
diff results/testBam2FastQReadName_2.fastq expected/testBam2FastQReadName_2.fastq
let "status |= $?"
diff results/testBam2FastQReadName.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoord 2> results/testBam2FastQCoord.log
let "status |= $?"
diff results/testBam2FastQCoord.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQCoord_1.fastq expected/testBam2FastQCoord_1.fastq
let "status |= $?"
diff results/testBam2FastQCoord_2.fastq expected/testBam2FastQCoord_2.fastq
let "status |= $?"
diff results/testBam2FastQCoord.log expected/testBam2FastQCoord.log
let "status |= $?"

##########################################
# Put the Read Name on the Plus Line
# Test clipping files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadNamePlus --rnPlus 2> results/testBam2FastQReadNamePlus.log
let "status |= $?"
diff results/testBam2FastQReadNamePlus.fastq expected/testBam2FastQReadNamePlus.fastq
let "status |= $?"
diff results/testBam2FastQReadNamePlus_1.fastq expected/testBam2FastQReadNamePlus_1.fastq
let "status |= $?"
diff results/testBam2FastQReadNamePlus_2.fastq expected/testBam2FastQReadNamePlus_2.fastq
let "status |= $?"
diff results/testBam2FastQReadNamePlus.log expected/testBam2FastQReadNamePlus.log
let "status |= $?"

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoordPlus --rnPlus 2> results/testBam2FastQCoordPlus.log
let "status |= $?"
diff results/testBam2FastQCoordPlus.fastq expected/testBam2FastQCoordPlus.fastq
let "status |= $?"
diff results/testBam2FastQCoordPlus_1.fastq expected/testBam2FastQCoordPlus_1.fastq
let "status |= $?"
diff results/testBam2FastQCoordPlus_2.fastq expected/testBam2FastQCoordPlus_2.fastq
let "status |= $?"
diff results/testBam2FastQCoordPlus.log expected/testBam2FastQCoordPlus.log
let "status |= $?"


##########################################
# Test specifying the Read Name Extension, not reverse complimenting, and
# specifically specifying some of the output files.
# Test clipping files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadNameNoCompBase --firstOut results/testBam2FastQReadNameNoCompFirst.fastq --unpairedOut results/testBam2FastQReadNameNoCompUnpaired.fastq --noReverseComp --firstRNExt _1  --secondRNExt _2 --rnPlus 2> results/testBam2FastQReadNameNoComp.log
let "status |= $?"
diff results/testBam2FastQReadNameNoCompUnpaired.fastq expected/testBam2FastQReadNameNoCompUnpaired.fastq
let "status |= $?"
diff results/testBam2FastQReadNameNoCompFirst.fastq expected/testBam2FastQReadNameNoCompFirst.fastq
let "status |= $?"
diff results/testBam2FastQReadNameNoCompBase_2.fastq expected/testBam2FastQReadNameNoCompBase_2.fastq
let "status |= $?"
diff results/testBam2FastQReadNameNoComp.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoordNoCompBase --secondOut results/testBam2FastQCoordNoCompSecond.fastq --noReverseComp --firstRNExt _1  --secondRNExt _2 2> results/testBam2FastQCoordNoComp.log
let "status |= $?"
diff results/testBam2FastQCoordNoCompBase.fastq expected/testBam2FastQCoordNoCompBase.fastq
let "status |= $?"
diff results/testBam2FastQCoordNoCompBase_1.fastq expected/testBam2FastQCoordNoCompBase_1.fastq
let "status |= $?"
diff results/testBam2FastQCoordNoCompSecond.fastq expected/testBam2FastQCoordNoCompSecond.fastq
let "status |= $?"
diff results/testBam2FastQCoordNoComp.log expected/testBam2FastQCoord.log
let "status |= $?"

##########################################
# Test BAM file inputs
# Test clipping files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testBam2FastQReadName.bam --outBase results/testBam2FastQBamReadName 2> results/testBam2FastQBamReadName.log
let "status |= $?"
diff results/testBam2FastQBamReadName.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQBamReadName_1.fastq expected/testBam2FastQReadName_1.fastq
let "status |= $?"
diff results/testBam2FastQBamReadName_2.fastq expected/testBam2FastQReadName_2.fastq
let "status |= $?"
diff results/testBam2FastQBamReadName.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.bam --outBase results/testBam2FastQBamCoord 2> results/testBam2FastQBamCoord.log
let "status |= $?"
diff results/testBam2FastQBamCoord.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQBamCoord_1.fastq expected/testBam2FastQCoord_1.fastq
let "status |= $?"
diff results/testBam2FastQBamCoord_2.fastq expected/testBam2FastQCoord_2.fastq
let "status |= $?"
diff results/testBam2FastQBamCoord.log expected/testBam2FastQCoord.log
let "status |= $?"

##########################################
# Test with files that specify sortorder
# Test clipping files sorted by read name.
../bin/bam bam2FastQ --in testFiles/testBam2FastQReadNameSO.sam --outBase results/testBam2FastQReadNameSO 2> results/testBam2FastQReadNameSO.log
let "status |= $?"
diff results/testBam2FastQReadNameSO.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQReadNameSO_1.fastq expected/testBam2FastQReadName_1.fastq
let "status |= $?"
diff results/testBam2FastQReadNameSO_2.fastq expected/testBam2FastQReadName_2.fastq
let "status |= $?"
diff results/testBam2FastQReadNameSO.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoordSO.sam --outBase results/testBam2FastQCoordSO 2> results/testBam2FastQCoordSO.log
let "status |= $?"
diff results/testBam2FastQCoordSO.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQCoordSO_1.fastq expected/testBam2FastQCoord_1.fastq
let "status |= $?"
diff results/testBam2FastQCoordSO_2.fastq expected/testBam2FastQCoord_2.fastq
let "status |= $?"
diff results/testBam2FastQCoordSO.log expected/testBam2FastQCoord.log
let "status |= $?"



if [ $status != 0 ]
then
  echo failed testBam2FastQ.sh
  exit 1
fi

