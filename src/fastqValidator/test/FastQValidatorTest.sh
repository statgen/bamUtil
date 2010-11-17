PATH_TO_EXE="../../bin"

ERROR=false

# Test running on a non-fastq file
$PATH_TO_EXE/fastQValidator --params --file ../FastQValidator.cpp > results/nonFastQFileResults.txt 2>&1


# Run on test file that tests all types of errors detected by the
# FastQValidator.  Specify to autodetect space type from the file
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --auto --printableErrors 100 --baseComposition > results/runResults.txt 2>&1
diff results/runResults.txt expectedResults/ExpectedResultsAutoDetect.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same file but do not check for unique sequence id.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --auto --printableErrors 100 --baseComposition --disableSeqIDCheck > results/runResultsDisableSeqID.txt 2>&1
diff results/runResultsDisableSeqID.txt expectedResults/ExpectedResultsDisableSeqID.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but specify Base Sequences.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 100 --baseComposition > results/runResultsBase.txt 2>&1
diff results/runResultsBase.txt expectedResults/ExpectedResultsBase.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but specify quit after -1 errors (do not quit until all read).
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 100 --baseComposition --maxErrors -1 > results/runResultsBaseNoQuit.txt 2>&1
diff results/runResultsBaseNoQuit.txt expectedResults/ExpectedResultsBase.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but specify quit after 0 errors (do not read the file).
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 100 --baseComposition --maxErrors 0 > results/runResultsBaseQuit0.txt 2>&1
diff results/runResultsBaseQuit0.txt expectedResults/ExpectedResults0Errors.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but specify quit after 10 errors, with only reporting the first 5 errors.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 5 --baseComposition --maxErrors 10 > results/runResultsBaseQuit10Err5.txt 2>&1
diff results/runResultsBaseQuit10Err5.txt expectedResults/ExpectedResults10Errors5Report.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but specify quit after 10 errors, reporting the first 100 errors, but only 10 will be reported because then it quit.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 100 --baseComposition --maxErrors 10 > results/runResultsBaseQuit10.txt 2>&1
diff results/runResultsBaseQuit10.txt expectedResults/ExpectedResults10Errors.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but specify Base Sequences, do not print baseComp.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --printableErrors 100 > results/runResultsNoBaseComp.txt 2>&1
diff results/runResultsNoBaseComp.txt expectedResults/ExpectedResultsNoBaseComp.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but ignore all errors, do not print baseComp.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --ignoreErrors > results/runResultsNoBaseCompNoErrors.txt 2>&1
diff results/runResultsNoBaseCompNoErrors.txt expectedResults/ExpectedResultsNoBaseCompNoErrors.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but do not print any messages, print baseComp.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --baseSpace --quiet --baseComposition > results/runResultsBaseCompNoMessages.txt 2>&1
diff results/runResultsBaseCompNoMessages.txt expectedResults/ExpectedResultsBaseCompNoMessages.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but auto detect as default.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --printableErrors 100 --baseComposition > results/runResultsBaseDefault.txt 2>&1
diff results/runResultsBaseDefault.txt expectedResults/ExpectedResultsAutoDetect.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but only accept Color Space Sequences
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --colorSpace --printableErrors 300 --baseComposition > results/runResultsColor.txt 2>&1
diff results/runResultsColor.txt expectedResults/ExpectedResultsColor.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on the same test file, but only accept Color Space Sequences limit
# number of printed errors to 5.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 10 --colorSpace --printableErrors 5 --baseComposition > results/runResultsColorLimitError.txt 2>&1
diff results/runResultsColorLimitError.txt expectedResults/ExpectedResultsColorLimitError.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

# Run on test file that tests all types of errors detected by the
# FastQValidator.  Accept Base Space as default.
# Modify to set minimum read length to 2.
$PATH_TO_EXE/fastQValidator --params --file testFile.txt --minReadLen 2 --printableErrors 100 --baseComposition > results/runResultsMinRead2.txt 2>&1
diff results/runResultsMinRead2.txt expectedResults/ExpectedResultsMinRead2.txt
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi

