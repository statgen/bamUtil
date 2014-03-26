# Sam To Sam, Reverting cigar & Qual
../bin/bam revert --in testFiles/testRevert.sam --out results/revertSam.sam --cigar --qual --noph 2> results/revertSam.sam.log && diff results/revertSam.sam expected/revertSam.sam && diff results/revertSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting cigar & Qual
../bin/bam revert --in testFiles/testRevert.sam --out results/revertSam.bam --cigar --qual --noph 2> results/revertSam.bam.log && diff results/revertSam.bam expected/revert.bam && diff results/revertSam.bam.log expected/revert.log \
&& \
# Sam To Sam, Reverting just cigar
../bin/bam revert --in testFiles/testRevert.sam --out results/revertCigarSam.sam --cigar --noph 2> results/revertCigarSam.sam.log && diff results/revertCigarSam.sam expected/revertCigarSam.sam && diff results/revertCigarSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting just cigar
../bin/bam revert --in testFiles/testRevert.sam --out results/revertCigarSam.bam --cigar --noph 2> results/revertCigarSam.bam.log && diff results/revertCigarSam.bam expected/revertCigar.bam && diff results/revertCigarSam.bam.log expected/revert.log \
&& \
# Sam To Sam, Reverting just qual
../bin/bam revert --in testFiles/testRevert.sam --out results/revertQualSam.sam --qual --noph 2> results/revertQualSam.sam.log && diff results/revertQualSam.sam expected/revertQualSam.sam && diff results/revertQualSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting just qual
../bin/bam revert --in testFiles/testRevert.sam --out results/revertQualSam.bam --qual --noph 2> results/revertQualSam.bam.log && diff results/revertQualSam.bam expected/revertQual.bam && diff results/revertQualSam.bam.log expected/revert.log \
&& \
# Sam To Sam, Reverting cigar & Qual without deleting them.
../bin/bam revert --in testFiles/testRevert.sam --out results/revertKeepTagsSam.sam --qual --cigar --keepTags --noph 2> results/revertKeepTagsSam.sam.log && diff results/revertKeepTagsSam.sam expected/revertKeepTagsSam.sam && diff results/revertKeepTagsSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting cigar & Qual without deleting them
../bin/bam revert --in testFiles/testRevert.sam --out results/revertKeepTagsSam.bam --qual --cigar --keepTags --noph 2> results/revertKeepTagsSam.bam.log && diff results/revertKeepTagsSam.bam expected/revertKeepTags.bam && diff results/revertKeepTagsSam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting cigar & Qual
../bin/bam revert --in testFiles/testRevert.bam --out results/revertBam.sam --cigar --qual --noph 2> results/revertBam.sam.log && diff results/revertBam.sam expected/revertBam.sam && diff results/revertBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting cigar & Qual
../bin/bam revert --in testFiles/testRevert.bam --out results/revertBam.bam --cigar --qual --noph 2> results/revertBam.bam.log && diff results/revertBam.bam expected/revert.bam && diff results/revertBam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting just cigar
../bin/bam revert --in testFiles/testRevert.bam --out results/revertCigarBam.sam --cigar --noph 2> results/revertCigarBam.sam.log && diff results/revertCigarBam.sam expected/revertCigarBam.sam && diff results/revertCigarBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting just cigar
../bin/bam revert --in testFiles/testRevert.bam --out results/revertCigarBam.bam --cigar --noph 2> results/revertCigarBam.bam.log && diff results/revertCigarBam.bam expected/revertCigar.bam && diff results/revertCigarBam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting just qual
../bin/bam revert --in testFiles/testRevert.bam --out results/revertQualBam.sam --qual --noph 2> results/revertQualBam.sam.log && diff results/revertQualBam.sam expected/revertQualBam.sam && diff results/revertQualBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting just qual
../bin/bam revert --in testFiles/testRevert.bam --out results/revertQualBam.bam --qual --noph 2> results/revertQualBam.bam.log && diff results/revertQualBam.bam expected/revertQual.bam && diff results/revertQualBam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting cigar & Qual without deleting them.
../bin/bam revert --in testFiles/testRevert.bam --out results/revertKeepTagsBam.sam --qual --cigar --keepTags --noph 2> results/revertKeepTagsBam.sam.log && diff results/revertKeepTagsBam.sam expected/revertKeepTagsBam.sam && diff results/revertKeepTagsBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting cigar & Qual without deleting them
../bin/bam revert --in testFiles/testRevert.bam --out results/revertKeepTagsBam.bam --qual --cigar --keepTags --noph 2> results/revertKeepTagsBam.bam.log && diff results/revertKeepTagsBam.bam expected/revertKeepTags.bam && diff results/revertKeepTagsBam.bam.log expected/revert.log \
&& \
\
# Sam To Sam, Reverting cigar & Qual, and deleting BQ
../bin/bam revert --in testFiles/testRevert.sam --out results/revertRmBQSam.sam --qual --cigar --rmBQ --noph 2> results/revertRmBQSam.sam.log && diff results/revertRmBQSam.sam expected/revertRmBQSam.sam && diff results/revertRmBQSam.sam.log expected/revert.log \
&& \
# Bam To Sam, Reverting cigar & Qual, and deleting BQ
../bin/bam revert --in testFiles/testRevert.bam --out results/revertRmBQBam.sam --qual --cigar --rmBQ --noph 2> results/revertRmBQBam.sam.log && diff results/revertRmBQBam.sam expected/revertRmBQBam.sam && diff results/revertRmBQBam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting cigar & Qual, and deleting BQ
../bin/bam revert --in testFiles/testRevert.sam --out results/revertRmBQSam.bam --qual --cigar --rmBQ --noph 2> results/revertRmBQSam.bam.log && diff results/revertRmBQSam.bam expected/revertRmBQSam.bam && diff results/revertRmBQSam.bam.log expected/revert.log \
&& \
# Bam To Bam, Reverting cigar & Qual, and deleting BQ
../bin/bam revert --in testFiles/testRevert.bam --out results/revertRmBQBam.bam --qual --cigar --rmBQ --noph 2> results/revertRmBQBam.bam.log && diff results/revertRmBQBam.bam expected/revertRmBQBam.bam && diff results/revertRmBQBam.bam.log expected/revert.log \
&& \
\
# Sam To Sam, Reverting cigar & Qual, and deleting MD:Z;AM:i
../bin/bam revert --in testFiles/testRevert.sam --out results/revertRmTagsSam.sam --qual --cigar --rmTags "MD:Z;AM:i" --noph 2> results/revertRmTagsSam.sam.log && diff results/revertRmTagsSam.sam expected/revertRmTagsSam.sam && diff results/revertRmTagsSam.sam.log expected/revert.log \
&& \
# Bam To Sam, Reverting cigar & Qual, and deleting MD:Z;AM:i;
../bin/bam revert --in testFiles/testRevert.bam --out results/revertRmTagsBam.sam --qual --cigar --rmTags "MD:Z;AM:i" --noph 2> results/revertRmTagsBam.sam.log && diff results/revertRmTagsBam.sam expected/revertRmTagsBam.sam && diff results/revertRmTagsBam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting cigar & Qual, and deleting MD:Z,AM:i;
../bin/bam revert --in testFiles/testRevert.sam --out results/revertRmTagsSam.bam --qual --cigar --rmTags "MD:Z,AM:i" --noph 2> results/revertRmTagsSam.bam.log && diff results/revertRmTagsSam.bam expected/revertRmTagsSam.bam && diff results/revertRmTagsSam.bam.log expected/revert.log \
&& \
# Bam To Bam, Reverting cigar & Qual, and deleting MD:Z;AM:i
../bin/bam revert --in testFiles/testRevert.bam --out results/revertRmTagsBam.bam --qual --cigar --rmTags "MD:Z;AM:i" --noph 2> results/revertRmTagsBam.bam.log && diff results/revertRmTagsBam.bam expected/revertRmTagsBam.bam && diff results/revertRmTagsBam.bam.log expected/revert.log \


if [ $? -ne 0 ]
then
  exit 1
fi

