#!/bin/bash

status=0;
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadName --noph 2> results/testBam2FastQReadName.log
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
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoord --noph 2> results/testBam2FastQCoord.log
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
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadNamePlus --rnPlus --noph 2> results/testBam2FastQReadNamePlus.log
let "status |= $?"
diff results/testBam2FastQReadNamePlus.fastq expected/testBam2FastQReadNamePlus.fastq
let "status |= $?"
diff results/testBam2FastQReadNamePlus_1.fastq expected/testBam2FastQReadNamePlus_1.fastq
let "status |= $?"
diff results/testBam2FastQReadNamePlus_2.fastq expected/testBam2FastQReadNamePlus_2.fastq
let "status |= $?"
diff results/testBam2FastQReadNamePlus.log expected/testBam2FastQReadNamePlus.log
let "status |= $?"

# Test converting files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoordPlus --rnPlus --noph 2> results/testBam2FastQCoordPlus.log
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
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadNameNoCompBase --firstOut results/testBam2FastQReadNameNoCompFirst.fastq --unpairedOut results/testBam2FastQReadNameNoCompUnpaired.fastq --noReverseComp --firstRNExt _1  --secondRNExt _2 --rnPlus --noph 2> results/testBam2FastQReadNameNoComp.log
let "status |= $?"
diff results/testBam2FastQReadNameNoCompUnpaired.fastq expected/testBam2FastQReadNameNoCompUnpaired.fastq
let "status |= $?"
diff results/testBam2FastQReadNameNoCompFirst.fastq expected/testBam2FastQReadNameNoCompFirst.fastq
let "status |= $?"
diff results/testBam2FastQReadNameNoCompBase_2.fastq expected/testBam2FastQReadNameNoCompBase_2.fastq
let "status |= $?"
diff results/testBam2FastQReadNameNoComp.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test converting files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoordNoCompBase --secondOut results/testBam2FastQCoordNoCompSecond.fastq --noReverseComp --firstRNExt _1  --secondRNExt _2 --noph 2> results/testBam2FastQCoordNoComp.log
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
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testBam2FastQReadName.bam --outBase results/testBam2FastQBamReadName --noph 2> results/testBam2FastQBamReadName.log
let "status |= $?"
diff results/testBam2FastQBamReadName.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQBamReadName_1.fastq expected/testBam2FastQReadName_1.fastq
let "status |= $?"
diff results/testBam2FastQBamReadName_2.fastq expected/testBam2FastQReadName_2.fastq
let "status |= $?"
diff results/testBam2FastQBamReadName.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test converting files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.bam --outBase results/testBam2FastQBamCoord --noph 2> results/testBam2FastQBamCoord.log
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
# Test with files that specify sortorder with reference for converting '='
# Test converting files sorted by read name.
../bin/bam bam2FastQ --in testFiles/testBam2FastQReadNameSO.sam --refFile testFilesLibBam/chr1_partial.fa --outBase results/testBam2FastQReadNameSO --noph 2> results/testBam2FastQReadNameSO.log
let "status |= $?"
diff results/testBam2FastQReadNameSO.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQReadNameSO_1.fastq expected/testBam2FastQReadNameSO_1.fastq
let "status |= $?"
diff results/testBam2FastQReadNameSO_2.fastq expected/testBam2FastQReadNameSO_2.fastq
let "status |= $?"
diff results/testBam2FastQReadNameSO.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test converting files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoordSO.sam --outBase results/testBam2FastQCoordSO --noph 2> results/testBam2FastQCoordSO.log
let "status |= $?"
diff results/testBam2FastQCoordSO.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQCoordSO_1.fastq expected/testBam2FastQCoord_1.fastq
let "status |= $?"
diff results/testBam2FastQCoordSO_2.fastq expected/testBam2FastQCoord_2.fastq
let "status |= $?"
diff results/testBam2FastQCoordSO.log expected/testBam2FastQCoord.log
let "status |= $?"


##########################################
# Test writing to just 2nd file to stdout.
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadName2ndStdout --secondOut - --noph 2> results/testBam2FastQReadName2ndStdout.log > results/testBam2FastQRN2ndStdout.fastq
let "status |= $?"
diff results/testBam2FastQReadName2ndStdout.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQReadName2ndStdout_1.fastq expected/testBam2FastQReadName_1.fastq
let "status |= $?"
diff results/testBam2FastQRN2ndStdout.fastq expected/testBam2FastQReadName_2.fastq
let "status |= $?"
diff results/testBam2FastQReadName2ndStdout.log expected/testBam2FastQReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoord2ndStdout --secondOut - --noph 2> results/testBam2FastQCoord2ndStdout.log > results/testBam2FastQCo2ndStdout.fastq
let "status |= $?"
diff results/testBam2FastQCoord2ndStdout.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQCoord2ndStdout_1.fastq expected/testBam2FastQCoord_1.fastq
let "status |= $?"
diff results/testBam2FastQCo2ndStdout.fastq expected/testBam2FastQCoord_2.fastq
let "status |= $?"
diff results/testBam2FastQCoord2ndStdout.log expected/testBam2FastQCoord.log
let "status |= $?"

##########################################
# Test merged output file for paired-end
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadNameMerge --merge --noph 2> results/testBam2FastQReadNameMerge.log
let "status |= $?"
diff results/testBam2FastQReadNameMerge.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQReadNameMerge_interleaved.fastq expected/testBam2FastQReadName_interleaved.fastq
let "status |= $?"
diff results/testBam2FastQReadNameMerge.log expected/testBam2FastQReadName.log
let "status |= $?"
if [ -e results/testBam2FastQReadNameMerge_1.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQReadNameMerge_2.fastq ]
then
  let "status = 1"
fi

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoordMerge --merge --noph 2> results/testBam2FastQCoordMerge.log
let "status |= $?"
diff results/testBam2FastQCoordMerge.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQCoordMerge_interleaved.fastq expected/testBam2FastQCoord_interleaved.fastq
let "status |= $?"
diff results/testBam2FastQCoordMerge.log expected/testBam2FastQCoord.log
let "status |= $?"
if [ -e results/testBam2FastQCoordMerge_1.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQCoordMerge_2.fastq ]
then
  let "status = 1"
fi


##########################################
# Test merged output file for paired-end by filenames rather than --merge
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --outBase results/testBam2FastQReadNameM --firstOut results/testBam2FastQRNM.fastq --secondOut results/testBam2FastQRNM.fastq --noph 2> results/testBam2FastQReadNameM.log
let "status |= $?"
diff results/testBam2FastQReadNameM.fastq expected/testBam2FastQReadName.fastq
let "status |= $?"
diff results/testBam2FastQRNM.fastq expected/testBam2FastQReadName_interleaved.fastq
let "status |= $?"
diff results/testBam2FastQReadNameM.log expected/testBam2FastQReadName.log
let "status |= $?"
if [ -e results/testBam2FastQReadNameM_1.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQReadNameM_2.fastq ]
then
  let "status = 1"
fi

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --outBase results/testBam2FastQCoordM --firstOut results/testBam2FastQCM.fastq --secondOut results/testBam2FastQCM.fastq --noph 2> results/testBam2FastQCoordM.log
let "status |= $?"
diff results/testBam2FastQCoordM.fastq expected/testBam2FastQCoord.fastq
let "status |= $?"
diff results/testBam2FastQCM.fastq expected/testBam2FastQCoord_interleaved.fastq
let "status |= $?"
diff results/testBam2FastQCoordM.log expected/testBam2FastQCoord.log
let "status |= $?"
if [ -e results/testBam2FastQCoordM_1.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQCoordM_2.fastq ]
then
  let "status = 1"
fi

##########################################
# Test merge all into same location
# Test converting files sorted by read name.
../bin/bam bam2FastQ --readName --in testFiles/testClipOverlapReadName.sam --unpair results/testBam2FastQReadNameAll1.fastq --firstOut results/testBam2FastQReadNameAll1.fastq --merge --noph 2> results/testBam2FastQReadNameAll1.log
let "status |= $?"
diff results/testBam2FastQReadNameAll1.fastq expected/testBam2FastQReadNameAll1.fastq
let "status |= $?"
diff results/testBam2FastQReadNameAll1.log expected/testBam2FastQReadName.log
let "status |= $?"
if [ -e results/testBam2FastQReadNameAll1_1.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQReadNameAll1_2.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQReadNameAll1_interleaved.fastq ]
then
  let "status = 1"
fi

# Test clipping files sorted by coordinate.
../bin/bam bam2FastQ --in testFiles/testBam2FastQCoord.sam --unpair results/testBam2FastQCoordAll1.fastq --firstOut results/testBam2FastQCoordAll1.fastq --merge --noph 2> results/testBam2FastQCoordAll1.log
let "status |= $?"
diff results/testBam2FastQCoordAll1.fastq expected/testBam2FastQCoordAll1.fastq
let "status |= $?"
diff results/testBam2FastQCoordAll1.log expected/testBam2FastQCoord.log
let "status |= $?"
if [ -e results/testBam2FastQCoordAll1_1.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQCoordAll1_2.fastq ]
then
  let "status = 1"
fi
if [ -e results/testBam2FastQCoordAll1_interleaved.fastq ]
then
  let "status = 1"
fi



if [ $status != 0 ]
then
  echo failed testBam2FastQ.sh
  exit 1
fi

