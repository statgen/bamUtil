ERROR=false

# squeeze sam to sam
../bin/bam squeeze --in testFiles/squeeze.sam --out results/squeezeSam.sam --noph 2> results/squeezeSam2Sam.log && \
diff results/squeezeSam.sam expected/squeeze.sam && diff results/squeezeSam2Sam.log expected/squeezeSam2Sam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to bam
../bin/bam squeeze --in testFiles/squeeze.bam --out results/squeezeBam.bam --noph 2> results/squeezeBam2Bam.log && \
diff results/squeezeBam.bam expected/squeeze.bam && diff results/squeezeBam2Bam.log expected/squeezeBam2Bam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to bam
../bin/bam squeeze --in testFiles/squeeze.sam --out results/squeezeSam.bam --noph 2> results/squeezeSam2Bam.log && \
diff results/squeezeSam.bam expected/squeeze.bam && diff results/squeezeSam2Bam.log expected/squeezeSam2Bam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to sam
../bin/bam squeeze --in testFiles/squeeze.bam --out results/squeezeBam.sam --noph 2> results/squeezeBam2Sam.log && \
diff results/squeezeBam.sam expected/squeezeBam.sam && diff results/squeezeBam2Sam.log expected/squeezeBam2Sam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to sam, keep OQ, keep Dups
../bin/bam squeeze --in testFiles/squeeze.sam --out results/squeezeKeep.sam --keepDups --keepOQ --refFile testFilesLibBam/chr1_partial.fa --noph 2> results/squeezeKeep.log && \
diff results/squeezeKeep.sam expected/squeezeKeep.sam && diff results/squeezeKeep.log expected/squeezeKeep.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using a low (default) quality value (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQual.sam --keepDups --keepOQ --binQualS 19,23,26,27 --noph 2> results/squeezeBinQual.log && \
diff results/squeezeBinQual.sam expected/squeezeBinQual.sam && diff results/squeezeBinQual.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using a high quality value (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQualHigh.sam --keepDups --keepOQ --binQualS 19,23,26,27 --binHigh --noph 2> results/squeezeBinQualHigh.log && \
diff results/squeezeBinQualHigh.sam expected/squeezeBinQualHigh.sam && diff results/squeezeBinQualHigh.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using the mid quality value (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQualMid.sam --keepDups --keepOQ --binQualS 19,23,26,27 --binMid --noph 2> results/squeezeBinQualMid.log && \
diff results/squeezeBinQualMid.sam expected/squeezeBinQualMid.sam && diff results/squeezeBinQualMid.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just binning qualities using the mid quality value, with range specified from a file (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeBinQualFileMid.sam --binMid --keepDups --keepOQ --binQualF testFiles/squeezeQualBin.txt --noph 2> results/squeezeBinQualFileMid.log && \
diff results/squeezeBinQualFileMid.sam expected/squeezeBinQualMid.sam && diff results/squeezeBinQualFileMid.log expected/squeezeBinQual.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sam to sam, just reducing read names (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeReadName.sam --readName results/squeezeReadNameMap.txt --keepDups --keepOQ --noph 2> results/squeezeReadName.log && \
diff results/squeezeReadName.sam expected/squeezeReadName.sam && diff results/squeezeReadName.log expected/squeezeReadName.log && diff results/squeezeReadNameMap.txt expected/squeezeReadNameMap.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to bam, just reducing read names (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testBam.bam --out results/squeezeReadName.bam --readName results/squeezeReadNameMapBam.txt --keepDups --keepOQ --noph 2> results/squeezeReadNameBam.log && \
diff results/squeezeReadName.bam expected/squeezeReadName.bam && diff results/squeezeReadNameBam.log expected/squeezeReadName.log && diff results/squeezeReadNameMapBam.txt expected/squeezeReadNameMap.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sorted by readname bam to sam, just reducing read names (keep OQ, keep dups).
../bin/bam squeeze --in testFiles/sortedReadName.sam --out results/squeezeReadNameSorted.sam --sreadName results/squeezeReadNameMapSamSorted.txt --keepDups --keepOQ --noph 2> results/squeezeReadNameSamSorted.log && \
diff results/squeezeReadNameSorted.sam expected/squeezeReadNameSorted.sam && diff results/squeezeReadNameSamSorted.log expected/squeezeReadName.log && diff results/squeezeReadNameMapSamSorted.txt expected/squeezeReadNameMapSorted.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze unsorted by read name bam to sam using the sorted method, just reducing read names (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testBam.bam --out results/squeezeReadNameUnsorted.sam --sreadName results/squeezeReadNameMapSamUnsorted.txt --keepDups --keepOQ --noph 2> results/squeezeReadNameSamUnsorted.log && \
diff results/squeezeReadNameUnsorted.sam expected/squeezeReadNameUnsorted.sam && diff results/squeezeReadNameSamUnsorted.log expected/squeezeReadName.log && diff results/squeezeReadNameMapSamUnsorted.txt expected/squeezeReadNameMapUnsorted.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze sorted by readname bam to bam, just reducing read names (keep OQ, keep dups).
../bin/bam squeeze --in testFiles/sortedReadName.bam --out results/squeezeReadNameSorted.bam --sreadName results/squeezeReadNameMapBamSorted.txt --keepDups --keepOQ --noph 2> results/squeezeReadNameBamSorted.log && \
diff results/squeezeReadNameSorted.bam expected/squeezeReadNameSorted.bam && diff results/squeezeReadNameBamSorted.log expected/squeezeReadName.log && diff results/squeezeReadNameMapBamSorted.txt expected/squeezeReadNameMapSorted.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze unsorted by read name bam to bam using the sorted method, just reducing read names (keep OQ, keep dups).
../bin/bam squeeze --in testFilesLibBam/testBam.bam --out results/squeezeReadNameUnsorted.bam --sreadName results/squeezeReadNameMapBamUnsorted.txt --keepDups --keepOQ --noph 2> results/squeezeReadNameBamUnsorted.log && \
diff results/squeezeReadNameUnsorted.bam expected/squeezeReadNameUnsorted.bam && diff results/squeezeReadNameBamUnsorted.log expected/squeezeReadName.log && diff results/squeezeReadNameMapBamUnsorted.txt expected/squeezeReadNameMapUnsorted.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to sam, keep dups, remove some tags.
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeTags.sam --keepDups --rmTags "XT:A;MD:Z;NM:i" --noph 2> results/squeezeTags.log && \
diff results/squeezeTags.sam expected/squeezeTags.sam && diff results/squeezeTags.log expected/squeezeTags.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to sam, keep dups, remove some tags.
../bin/bam squeeze --in testFilesLibBam/testSam.sam --out results/squeezeTags1.sam --keepDups --rmTags "XT:A,MD:Z,NM:i" --noph 2> results/squeezeTags1.log && \
diff results/squeezeTags1.sam expected/squeezeTags.sam && diff results/squeezeTags1.log expected/squeezeTags.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
    echo "Fail testSqueeze.sh"
  exit 1
fi

