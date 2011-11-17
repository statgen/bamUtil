#!/bin/bash

status=0;
# Test clipping files sorted by read name.
../bin/bam clipOverlap --readName --in testFiles/testClipOverlapReadName.sam --out results/testClipOverlapReadName.sam --storeOrig XC 2> results/testClipOverlapReadName.log
let "status |= $?"
diff results/testClipOverlapReadName.sam expected/testClipOverlapReadName.sam
let "status |= $?"
diff results/testClipOverlapReadName.log expected/testClipOverlapReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoord.sam --storeOrig XC 2> results/testClipOverlapCoord.log
let "status |= $?"
diff results/testClipOverlapCoord.sam expected/testClipOverlapCoord.sam
let "status |= $?"
diff results/testClipOverlapCoord.log expected/testClipOverlapCoord.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool without default clipping
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3.sam --storeOrig XC --poolSize 3 --poolSkipClip 2> results/testClipOverlapCoordPool3.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3.sam expected/testClipOverlapCoordPool3.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3.log expected/testClipOverlapCoordPool3.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool with default clipping
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3Clip.sam --storeOrig XC --poolSize 3 2> results/testClipOverlapCoordPool3Clip.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3Clip.sam expected/testClipOverlapCoordPool3Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3Clip.log expected/testClipOverlapCoordPool3Clip.log
let "status |= $?"

# Test clipping files sorted by coordinate with no pool with default clipping
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool0Clip.sam --storeOrig XC --poolSize 0 2> results/testClipOverlapCoordPool0Clip.log
if [ $? != 8 ]
then
    status=1
    echo did not get expected return value for no pool.
fi
diff results/testClipOverlapCoordPool0Clip.sam expected/testClipOverlapCoordPool0Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool0Clip.log expected/testClipOverlapCoordPool0Clip.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed testClipOverlap.sh
  exit 1
fi

