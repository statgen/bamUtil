#!/bin/bash

status=0;
../bin/bam dedup --in testFiles/testDedup.sam --out results/testDedup.sam 2> results/testDedup.txt
let "status |= $?"
diff results/testDedup.txt expected/testDedup.txt
let "status |= $?"
diff results/testDedup.sam expected/testDedup.sam
let "status |= $?"


if [ $status != 0 ]
then
  echo failed testDedup.sh
  exit 1
fi

