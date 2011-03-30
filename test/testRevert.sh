# Sam To Sam, Reverting cigar & Qual
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertSam.sam --cigar --qual 2> results/revertSam.sam.log && diff results/revertSam.sam expected/revertSam.sam && diff results/revertSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting cigar & Qual
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertSam.bam --cigar --qual 2> results/revertSam.bam.log && diff results/revertSam.bam expected/revertSam.bam && diff results/revertSam.bam.log expected/revert.log \
&& \
# Sam To Sam, Reverting just cigar
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertCigarSam.sam --cigar 2> results/revertCigarSam.sam.log && diff results/revertCigarSam.sam expected/revertCigarSam.sam && diff results/revertCigarSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting just cigar
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertCigarSam.bam --cigar 2> results/revertCigarSam.bam.log && diff results/revertCigarSam.bam expected/revertCigarSam.bam && diff results/revertCigarSam.bam.log expected/revert.log \
&& \
# Sam To Sam, Reverting just qual
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertQualSam.sam --qual 2> results/revertQualSam.sam.log && diff results/revertQualSam.sam expected/revertQualSam.sam && diff results/revertQualSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting just qual
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertQualSam.bam --qual 2> results/revertQualSam.bam.log && diff results/revertQualSam.bam expected/revertQualSam.bam && diff results/revertQualSam.bam.log expected/revert.log \
&& \
# Sam To Sam, Reverting cigar & Qual without deleting them.
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertKeepTagsSam.sam --qual --cigar --keepTags 2> results/revertKeepTagsSam.sam.log && diff results/revertKeepTagsSam.sam expected/revertKeepTagsSam.sam && diff results/revertKeepTagsSam.sam.log expected/revert.log \
&& \
# Sam To Bam, Reverting cigar & Qual without deleting them
../../bin/bam revert --in testFiles/testRevert.sam --out results/revertKeepTagsSam.bam --qual --cigar --keepTags 2> results/revertKeepTagsSam.bam.log && diff results/revertKeepTagsSam.bam expected/revertKeepTagsSam.bam && diff results/revertKeepTagsSam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting cigar & Qual
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertBam.sam --cigar --qual 2> results/revertBam.sam.log && diff results/revertBam.sam expected/revertBam.sam && diff results/revertBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting cigar & Qual
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertBam.bam --cigar --qual 2> results/revertBam.bam.log && diff results/revertBam.bam expected/revertBam.bam && diff results/revertBam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting just cigar
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertCigarBam.sam --cigar 2> results/revertCigarBam.sam.log && diff results/revertCigarBam.sam expected/revertCigarBam.sam && diff results/revertCigarBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting just cigar
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertCigarBam.bam --cigar 2> results/revertCigarBam.bam.log && diff results/revertCigarBam.bam expected/revertCigarBam.bam && diff results/revertCigarBam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting just qual
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertQualBam.sam --qual 2> results/revertQualBam.sam.log && diff results/revertQualBam.sam expected/revertQualBam.sam && diff results/revertQualBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting just qual
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertQualBam.bam --qual 2> results/revertQualBam.bam.log && diff results/revertQualBam.bam expected/revertQualBam.bam && diff results/revertQualBam.bam.log expected/revert.log \
&& \
# Bam To Sam, Reverting cigar & Qual without deleting them.
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertKeepTagsBam.sam --qual --cigar --keepTags 2> results/revertKeepTagsBam.sam.log && diff results/revertKeepTagsBam.sam expected/revertKeepTagsBam.sam && diff results/revertKeepTagsBam.sam.log expected/revert.log \
&& \
# Bam To Bam, Reverting cigar & Qual without deleting them
../../bin/bam revert --in testFiles/testRevert.bam --out results/revertKeepTagsBam.bam --qual --cigar --keepTags 2> results/revertKeepTagsBam.bam.log && diff results/revertKeepTagsBam.bam expected/revertKeepTagsBam.bam && diff results/revertKeepTagsBam.bam.log expected/revert.log 



if [ $? -ne 0 ]
then
  exit 1
fi

