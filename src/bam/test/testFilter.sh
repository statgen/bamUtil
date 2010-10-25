ERROR=false

../../bin/bam filter --ref testFiles/chr1_partial.fa --in testFiles/testFilter.sam --mis .49 --qu 30  > results/filter.sam ; diff results/filter.sam expected/filter.sam

if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi

