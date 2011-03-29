ERROR=false

# squeeze sam to sam
../../bin/squeezeBam testFiles/squeeze.sam results/squeezeSam.sam 2> results/sam2sam.log
diff results/squeezeSam.sam expected/squeeze.sam && diff results/sam2sam.log expected/sam2sam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to bam
../../bin/squeezeBam testFiles/squeeze.bam results/squeezeBam.bam 2> results/bam2bam.log
diff results/squeezeBam.bam expected/squeeze.bam && diff results/bam2bam.log expected/bam2bam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


# squeeze sam to bam
../../bin/squeezeBam testFiles/squeeze.sam results/squeezeSam.bam 2> results/sam2bam.log
diff results/squeezeSam.bam expected/squeeze.bam && diff results/sam2bam.log expected/sam2bam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# squeeze bam to sam
../../bin/squeezeBam testFiles/squeeze.bam results/squeezeBam.sam 2> results/bam2sam.log
diff results/squeezeBam.sam expected/squeezeBam.sam && diff results/bam2sam.log expected/bam2sam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  exit 1
fi

