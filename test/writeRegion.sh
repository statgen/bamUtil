ERROR=false

../bin/bam  writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionNeg1.bam --refID -1 2> results/regionNeg1.txt \
&& diff results/regionNeg1.bam expected/regionNeg1.bam && diff results/regionNeg1.txt expected/regionNeg1.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead.sam --readName "1:1011:F:255+17M15D20M" 2> results/regionRead.txt \
&& diff results/regionRead.sam expected/regionRead.sam && diff results/regionRead.txt expected/regionRead.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead2.sam --readName "1:1011:F:255+17M15D20M" --refName 1 --start 1010 --end 1011 2> results/regionRead2.txt \
&& diff results/regionRead2.sam expected/regionRead2.sam && diff results/regionRead2.txt expected/regionRead2.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead3.sam --readName "1:1011:F:255+17M15D20" --refName 1 --start 1010 --end 1011 2> results/regionRead3.txt \
&& diff results/regionRead3.sam expected/regionReadNone.sam && diff results/regionRead3.txt expected/regionRead3.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead4.sam --refName 1 --start 1010 --end 1011 2> results/regionRead4.txt \
&& diff results/regionRead4.sam expected/regionRead2.sam && diff results/regionRead4.txt expected/regionRead4.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead5.sam --refName 1 --start 1010 --end 1011 --withinReg 2> results/regionRead5.txt \
&& diff results/regionRead5.sam expected/regionReadNone.sam && diff results/regionRead5.txt expected/regionRead5.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead6.sam --refName 1 --start 1014 --end 1015 2> results/regionRead6.txt \
&& diff results/regionRead6.sam expected/regionRead2.sam && diff results/regionRead6.txt expected/regionRead6.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead7.sam --refName 1 --start 1014 --end 1015 --withinReg 2> results/regionRead7.txt \
&& diff results/regionRead7.sam expected/regionReadNone.sam && diff results/regionRead7.txt expected/regionRead7.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead8.sam --bed testFiles/bedFile.bed --withinReg 2> results/regionRead8.txt \
&& diff results/regionRead8.sam expected/regionRead8.sam && diff results/regionRead8.txt expected/regionRead8.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead9.sam --bed testFiles/bedFile.bed 2> results/regionRead9.txt \
&& diff results/regionRead9.sam expected/regionRead9.sam && diff results/regionRead9.txt expected/regionRead9.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead10.sam --bed testFiles/bedFile2.bed --withinReg 2> results/regionRead10.txt \
&& diff results/regionRead10.sam expected/regionRead8.sam && diff results/regionRead10.txt expected/regionRead10.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionRead11.sam --bed testFiles/bedFile2.bed 2> results/regionRead11.txt \
&& diff results/regionRead11.sam expected/regionRead9.sam && diff results/regionRead11.txt expected/regionRead11.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/testShift.sam --out results/regionShift.sam --lshift 2> results/regionShift.txt \
&& diff results/regionShift.sam expected/regionShift.sam && diff results/regionShift.txt expected/regionShift.txt\
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionReadRnFile.sam --rnFile testFiles/rn.txt 2> results/regionReadRnFile.txt \
&& diff results/regionReadRnFile.sam expected/regionRead.sam && diff results/regionReadRnFile.txt expected/regionReadRnFile.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionReadRnFile1.sam --rnFile testFiles/rn1.txt 2> results/regionReadRnFile1.txt \
&& diff results/regionReadRnFile1.sam expected/regionRead.sam && diff results/regionReadRnFile1.txt expected/regionReadRnFile1.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionReadRnFile2.sam --rnFile testFiles/rn2.txt 2> results/regionReadRnFile2.txt \
&& diff results/regionReadRnFile2.sam expected/regionReadRnFile2.sam && diff results/regionReadRnFile2.txt expected/regionReadRnFile2.txt \
&& \
../bin/bam writeRegion --noph --in testFilesLibBam/sortedBam.bam --out results/regionReadRnFile3.sam --rnFile testFiles/rn3.txt 2> results/regionReadRnFile3.txt \
&& diff results/regionReadRnFile3.sam expected/regionReadRnFile2.sam && diff results/regionReadRnFile3.txt expected/regionReadRnFile3.txt \

if [ $? -ne 0 ]
then
  exit 1
fi

