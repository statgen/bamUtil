ERROR=false

../../bin/polishBam  --in testFiles/sortedSam.sam --out results/updatedSam.sam -l results/updated.log --checkSQ --fasta testFiles/testFasta.fa --AS my37 --UR testFasta.fa --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none"

../../bin/polishBam  -i testFiles/sortedSam.sam -o results/updatedSam1.sam -l results/updated1.log --checkSQ -f testFiles/testFasta.fa --AS my37 --UR testFasta.fa --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none"

diff results/updatedSam.sam expected/updatedSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/updatedSam1.sam expected/updatedSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi


if($ERROR == true)
then
  exit 1
fi
