ERROR=false

../bin/bam polishBam  --in testFiles/sortedSam.sam --out results/polishSam.sam --log results/polishSam.log --checkSQ --fasta testFiles/testFasta.fa --AS my37 --UR testFasta.fa --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none" --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

../bin/bam polishBam  -i testFiles/sortedSam.sam -o results/polishSam1.sam -l results/polishSam1.log --checkSQ -f testFiles/testFasta.fa --AS my37 --UR testFasta.fa --RG "@RG	ID:UM0037:1	SM:Sample2	LB:lb2	PU:mypu	CN:UMCORE	DT:2010-11-01	PL:ILLUMINA" --PG "@PG	ID:polish	VN:0.0.1" --SP new --HD "@HD	VN:1.0	SO:coordinate	GO:none" --noph
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/polishSam.sam expected/polishSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/polishSam.log expected/polishSam.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/polishSam1.sam expected/polishSam.sam
if [ $? -ne 0 ]
then
    ERROR=true
fi

diff results/polishSam1.log expected/polishSam1.log
if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi
