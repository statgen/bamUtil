ERROR=false

# Convert SAM to SAM
../bin/bam convert --in testFiles/testTags.sam --out results/testTagsSam2Sam.sam --noph > results/testTagsSam2Sam.log 2> results/testTagsSam2Sam.err 
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/testTagsSam2Sam.sam expected/testTags.sam && \
diff results/testTagsSam2Sam.log expected/empty.txt && \
diff results/testTagsSam2Sam.err expected/testTags.err
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Convert SAM to BAM
../bin/bam convert --in testFiles/testTags.sam --out results/testTagsSam2Bam.bam --noph > results/testTagsSam2Bam.log 2> results/testTagsSam2Bam.err 
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/testTagsSam2Bam.bam expected/testTags.bam && \
diff results/testTagsSam2Bam.log expected/empty.txt && \
diff results/testTagsSam2Bam.err expected/testTags.err
if [ $? -ne 0 ]
then
    ERROR=true
fi



# Convert BAM to SAM
../bin/bam convert --in testFiles/testTags.bam --out results/testTagsBam2Sam.sam --noph > results/testTagsBam2Sam.log 2> results/testTagsBam2Sam.err 
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/testTagsBam2Sam.sam expected/testTags.sam && \
diff results/testTagsBam2Sam.log expected/empty.txt && \
diff results/testTagsBam2Sam.err expected/testTags.err
if [ $? -ne 0 ]
then
    ERROR=true
fi




# Convert BAM to BAM - doesn't remove duplicates since the tags aren't parsed.
../bin/bam convert --in testFiles/testTags.bam --out results/testTagsBam2Bam.bam --noph > results/testTagsBam2Bam.log 2> results/testTagsBam2Bam.err 
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/testTagsBam2Bam.bam expected/testTagsBam2Bam.bam && \
diff results/testTagsBam2Bam.log expected/empty.txt && \
diff results/testTagsBam2Bam.err expected/testTagsBam2Bam.err
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  echo "Failed testTags.sh"
  exit 1
fi
