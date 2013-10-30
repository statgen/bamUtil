ERROR=false

../bin/bam filter --ref testFilesLibBam/chr1_partial.fa --in testFiles/testFilter.sam --mis .49 --qu 30 --noph > results/filter.sam 2> results/filter.log && diff results/filter.sam expected/filter.sam && diff results/filter.log expected/filter.log

if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi

