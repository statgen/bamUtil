ERROR=false

../bin/bam validate --params --in testFiles/testInvalid.sam  --v 2> results/validateInvalid.txt
diff results/validateInvalid.txt expected/invalid.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam validate --params --in testFiles/testInvalid1.sam  --v 2> results/validateInvalid1.txt
diff results/validateInvalid1.txt expected/invalid1.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting bam to sam 
../bin/bam convert --params --in testFilesLibBam/testBam.bam --out results/convertBam.sam 2> results/convertBam.log && diff results/convertBam.sam expected/convertBam.sam && diff results/convertBam.log expected/convertBam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting bam to sam, writing to stdout
../bin/bam convert --params --in testFilesLibBam/testBam.bam --out - > results/convertBamStdout.sam 2> results/convertBamStdout.log && diff results/convertBamStdout.sam expected/convertBam.sam && diff results/convertBamStdout.log expected/convertBamStdout.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting sam to bam 
../bin/bam convert --params --in testFilesLibBam/testSam.sam --out results/convertSam.bam 2> results/convertSam.log && diff results/convertSam.bam expected/convertSam.bam && diff results/convertSam.log expected/convertSam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting sam to bam, writing to stdout
../bin/bam convert --params --in testFilesLibBam/testSam.sam --out -.bam > results/convertSamStdout.bam 2> results/convertSamStdout.log && diff results/convertSamStdout.bam expected/convertSam.bam && diff results/convertSamStdout.log expected/convertSamStdout.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting sam to bam, writing to stdout
../bin/bam convert --params --in testFilesLibBam/testSam.sam --out -.bam > results/convertSamStdout.bam 2> results/convertSamStdout.log && diff results/convertSamStdout.bam expected/convertSam.bam && diff results/convertSamStdout.log expected/convertSamStdout.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting sam to bam and back to sam via stdout, pipe, and stdin reading bam via stdin
../bin/bam convert --params --in testFilesLibBam/testSam.sam --out -.bam  2> results/convertSamStdoutBamPipe.log | ../bin/bam convert --params --in -.bam --out results/convertSamStdoutBamPipeSam.sam 2> results/convertSamStdoutBamPipeSam.log && diff results/convertSamStdoutBamPipeSam.sam expected/convertBam.sam && diff results/convertSamStdoutBamPipeSam.log expected/convertSamStdoutBamPipeSam.log && diff results/convertSamStdoutBamPipe.log expected/convertSamStdoutBamPipe.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Test converting bam to sam and back to bam via stdout, pipe, and stdin reading sam via stdin
../bin/bam convert --params --in testFilesLibBam/testBam.bam --out -.sam  2> results/convertBamStdoutSamPipe.log | ../bin/bam convert --params --in -.sam --out results/convertBamStdoutSamPipeBam.bam 2> results/convertBamStdoutSamPipeBam.log && diff results/convertBamStdoutSamPipeBam.bam expected/convertSam.bam && diff results/convertBamStdoutSamPipeBam.log expected/convertBamStdoutSamPipeBam.log && diff results/convertBamStdoutSamPipe.log expected/convertBamStdoutSamPipe.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

#TODO test reading & writing ubam over stdin/stdout.


if($ERROR == true)
then
  exit 1
fi

