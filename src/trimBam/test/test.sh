ERROR=false

../../bin/trimBam testFiles/testSam.sam results/trimSam.sam 2 2> results/test.log
diff results/trimSam.sam expected/trimSam.sam && diff results/test.log expected/test.log
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  exit 1
fi

