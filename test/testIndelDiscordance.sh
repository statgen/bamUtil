#!/bin/bash

status=0;
# Test indelDiscordance
../bin/bam indelDiscordance --in testFiles/testIndelDiscordance.bam --refFile testFilesLibBam/chr1_partial.fa --chrom 1 --start 0 --end 100000 2> results/discordance.txt
let "status |= $?"
diff results/discordance.txt expected/discordance.txt
let "status |= $?"


if [ $status != 0 ]
then
  echo failed testIndelDiscordance.sh
  exit 1
fi

