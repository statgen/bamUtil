ERROR=false

../../bin/trimBam testFiles/testSam.sam results/trimSam.sam 2
diff results/trimSam.sam expected/trimSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  exit 1
fi

