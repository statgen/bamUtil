ERROR=false

# squeeze sam to sam
../bin/bam squeeze --in testFiles/squeeze.sam --out results/squeezeSam.sam 2> results/squeezeSam2Sam.log && \
diff results/squeezeSam.sam expected/squeeze.sam && diff results/squeezeSam2Sam.log expected/squeezeSam2Sam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to bam
../bin/bam squeeze --in testFiles/squeeze.bam --out results/squeezeBam.bam 2> results/squeezeBam2Bam.log && \
diff results/squeezeBam.bam expected/squeeze.bam && diff results/squeezeBam2Bam.log expected/squeezeBam2Bam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to bam
../bin/bam squeeze --in testFiles/squeeze.sam --out results/squeezeSam.bam 2> results/squeezeSam2Bam.log && \
diff results/squeezeSam.bam expected/squeeze.bam && diff results/squeezeSam2Bam.log expected/squeezeSam2Bam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to sam
../bin/bam squeeze --in testFiles/squeeze.bam --out results/squeezeBam.sam 2> results/squeezeBam2Sam.log && \
diff results/squeezeBam.sam expected/squeezeBam.sam && diff results/squeezeBam2Sam.log expected/squeezeBam2Sam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to sam, keep OQ, keep Dups
../bin/bam squeeze --in testFiles/squeeze.sam --out results/squeezeKeep.sam --keepDups --keepOQ --refFile testFilesLibBam/chr1_partial.fa 2> results/squeezeKeep.log && \
diff results/squeezeKeep.sam expected/squeezeKeep.sam && diff results/squeezeKeep.log expected/squeezeKeep.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using a low (default) quality value (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQual.sam --keepDups --keepOQ --binQualS 19,23,26,27 2> results/squeezeBinQual.log && \
diff results/squeezeBinQual.sam expected/squeezeBinQual.sam && diff results/squeezeBinQual.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using a high quality value (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQualHigh.sam --keepDups --keepOQ --binQualS 19,23,26,27 --binHigh 2> results/squeezeBinQualHigh.log && \
diff results/squeezeBinQualHigh.sam expected/squeezeBinQualHigh.sam && diff results/squeezeBinQualHigh.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using the mid quality value (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQualMid.sam --keepDups --keepOQ --binQualS 19,23,26,27 --binMid 2> results/squeezeBinQualMid.log && \
diff results/squeezeBinQualMid.sam expected/squeezeBinQualMid.sam && diff results/squeezeBinQualMid.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using the mid quality value, with range specified from a file (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQualFileMid.sam --binMid --keepDups --keepOQ --binQualF testFiles/squeezeQualBin.txt 2> results/squeezeBinQualFileMid.log && \
diff results/squeezeBinQualFileMid.sam expected/squeezeBinQualMid.sam && diff results/squeezeBinQualFileMid.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi

