#!/bin/bash

status=0;
../bin/bam gapInfo --in testFiles/testClipOverlapCoord.sam --out results/gapInfo.txt 2> results/gapInfo.log
let "status |= $?"
diff results/gapInfo.txt expected/gapInfo.txt
let "status |= $?"
diff results/gapInfo.log expected/empty.log
let "status |= $?"


if [ $status != 0 ]
then
  echo failed testClipOverlap.sh
  exit 1
fi

