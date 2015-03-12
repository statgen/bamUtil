#!/bin/bash

status=0;
# Test clipping files sorted by read name.
../bin/bam clipOverlap --readName --in testFiles/testClipOverlapReadName.sam --out results/testClipOverlapReadName.sam --storeOrig XC --noph 2> results/testClipOverlapReadName.log
let "status |= $?"
diff results/testClipOverlapReadName.sam expected/testClipOverlapReadName.sam
let "status |= $?"
diff results/testClipOverlapReadName.log expected/testClipOverlapReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoord.sam --storeOrig XC --noph 2> results/testClipOverlapCoord.log
let "status |= $?"
diff results/testClipOverlapCoord.sam expected/testClipOverlapCoord.sam
let "status |= $?"
diff results/testClipOverlapCoord.log expected/testClipOverlapCoord.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordUnmap.sam --storeOrig XC --noph --unmap 2> results/testClipOverlapCoordUnmap.log
let "status |= $?"
diff results/testClipOverlapCoordUnmap.sam expected/testClipOverlapCoordUnmap.sam
let "status |= $?"
diff results/testClipOverlapCoordUnmap.log expected/testClipOverlapCoord.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool without default clipping
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3.sam --storeOrig XC --poolSize 3 --poolSkipClip --noph 2> results/testClipOverlapCoordPool3.log
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
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3Clip.sam --storeOrig XC --poolSize 3 --noph 2> results/testClipOverlapCoordPool3Clip.log
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
../bin/bam clipOverlap --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool0Clip.sam --storeOrig XC --poolSize 0 --noph 2> results/testClipOverlapCoordPool0Clip.log
if [ $? != 8 ]
then
    status=1
    echo did not get expected return value for no pool.
fi
diff results/testClipOverlapCoordPool0Clip.sam expected/testClipOverlapCoordPool0Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool0Clip.log expected/testClipOverlapCoordPool0Clip.log
let "status |= $?"

# Test clipping files sorted by read name.
../bin/bam clipOverlap --stats --readName --in testFiles/testClipOverlapReadName.sam --out results/testClipOverlapReadNameStats.sam --storeOrig XC --noph 2> results/testClipOverlapReadNameStats.log
let "status |= $?"
diff results/testClipOverlapReadNameStats.sam expected/testClipOverlapReadName.sam
let "status |= $?"
diff results/testClipOverlapReadNameStats.log expected/testClipOverlapReadNameStats.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --stats --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordStats.sam --storeOrig XC --noph 2> results/testClipOverlapCoordStats.log
let "status |= $?"
diff results/testClipOverlapCoordStats.sam expected/testClipOverlapCoord.sam
let "status |= $?"
diff results/testClipOverlapCoordStats.log expected/testClipOverlapCoordStats.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool without default clipping
../bin/bam clipOverlap --stats --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3Stats.sam --storeOrig XC --poolSize 3 --poolSkipClip --noph 2> results/testClipOverlapCoordPool3Stats.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3Stats.sam expected/testClipOverlapCoordPool3.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3Stats.log expected/testClipOverlapCoordPool3Stats.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool with default clipping
../bin/bam clipOverlap --stats --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3ClipStats.sam --storeOrig XC --poolSize 3 --noph 2> results/testClipOverlapCoordPool3ClipStats.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3ClipStats.sam expected/testClipOverlapCoordPool3Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3ClipStats.log expected/testClipOverlapCoordPool3ClipStats.log
let "status |= $?"

# Test clipping files sorted by coordinate with no pool with default clipping
../bin/bam clipOverlap --stats --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool0ClipStats.sam --storeOrig XC --poolSize 0 --noph 2> results/testClipOverlapCoordPool0ClipStats.log
if [ $? != 8 ]
then
    status=1
    echo did not get expected return value for no pool.
fi
diff results/testClipOverlapCoordPool0ClipStats.sam expected/testClipOverlapCoordPool0Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool0ClipStats.log expected/testClipOverlapCoordPool0ClipStats.log
let "status |= $?"


#########################################
# Test with only keeping the clipped read pairs.
# Test clipping files sorted by read name.
../bin/bam clipOverlap --readName --clipsOnly --in testFiles/testClipOverlapReadName.sam --out results/testClipOverlapReadNameClipsOnly.sam --storeOrig XC --noph 2> results/testClipOverlapReadNameClipsOnly.log
let "status |= $?"
diff results/testClipOverlapReadNameClipsOnly.sam expected/testClipOverlapReadNameClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapReadNameClipsOnly.log expected/testClipOverlapReadName.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordClipsOnly.sam --storeOrig XC --noph 2> results/testClipOverlapCoordClipsOnly.log
let "status |= $?"
diff results/testClipOverlapCoordClipsOnly.sam expected/testClipOverlapCoordClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapCoordClipsOnly.log expected/testClipOverlapCoord.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool without default clipping
../bin/bam clipOverlap --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3ClipsOnly.sam --storeOrig XC --poolSize 3 --poolSkipClip --noph 2> results/testClipOverlapCoordPool3ClipsOnly.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3ClipsOnly.sam expected/testClipOverlapCoordPool3ClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3ClipsOnly.log expected/testClipOverlapCoordPool3ClipsOnly.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool with default clipping
../bin/bam clipOverlap --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3ClipClipsOnly.sam --storeOrig XC --poolSize 3 --noph 2> results/testClipOverlapCoordPool3ClipClipsOnly.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3ClipClipsOnly.sam expected/testClipOverlapCoordPool3ClipClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3ClipClipsOnly.log expected/testClipOverlapCoordPool3ClipClipsOnly.log
let "status |= $?"

# Test clipping files sorted by coordinate with no pool with default clipping
../bin/bam clipOverlap --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool0ClipClipsOnly.sam --storeOrig XC --poolSize 0 --noph 2> results/testClipOverlapCoordPool0ClipClipsOnly.log
if [ $? != 8 ]
then
    status=1
    echo did not get expected return value for no pool.
fi
diff results/testClipOverlapCoordPool0ClipClipsOnly.sam expected/testClipOverlapCoordPool0Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool0ClipClipsOnly.log expected/testClipOverlapCoordPool0Clip.log
let "status |= $?"

#### Add Stats Option
# Test clipping files sorted by read name.
../bin/bam clipOverlap --stats --readName --clipsOnly --in testFiles/testClipOverlapReadName.sam --out results/testClipOverlapReadNameStatsClipsOnly.sam --storeOrig XC --noph 2> results/testClipOverlapReadNameStatsClipsOnly.log
let "status |= $?"
diff results/testClipOverlapReadNameStatsClipsOnly.sam expected/testClipOverlapReadNameClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapReadNameStatsClipsOnly.log expected/testClipOverlapReadNameStats.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --stats --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordStatsClipsOnly.sam --storeOrig XC --noph 2> results/testClipOverlapCoordStatsClipsOnly.log
let "status |= $?"
diff results/testClipOverlapCoordStatsClipsOnly.sam expected/testClipOverlapCoordClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapCoordStatsClipsOnly.log expected/testClipOverlapCoordStats.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool without default clipping
../bin/bam clipOverlap --stats --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3StatsClipsOnly.sam --storeOrig XC --poolSize 3 --poolSkipClip --noph 2> results/testClipOverlapCoordPool3StatsClipsOnly.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3StatsClipsOnly.sam expected/testClipOverlapCoordPool3ClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3StatsClipsOnly.log expected/testClipOverlapCoordPool3StatsClipsOnly.log
let "status |= $?"

# Test clipping files sorted by coordinate with small pool with default clipping
../bin/bam clipOverlap --stats --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool3ClipStatsClipsOnly.sam --storeOrig XC --poolSize 3 --noph 2> results/testClipOverlapCoordPool3ClipStatsClipsOnly.log
if [ $? != 2 ]
then
    status=1
    echo did not get expected return value for a small pool.
fi
diff results/testClipOverlapCoordPool3ClipStatsClipsOnly.sam expected/testClipOverlapCoordPool3ClipClipsOnly.sam
let "status |= $?"
diff results/testClipOverlapCoordPool3ClipStatsClipsOnly.log expected/testClipOverlapCoordPool3ClipStatsClipsOnly.log
let "status |= $?"

# Test clipping files sorted by coordinate with no pool with default clipping
../bin/bam clipOverlap --stats --clipsOnly --in testFiles/testClipOverlapCoord.sam --out results/testClipOverlapCoordPool0ClipStatsClipsOnly.sam --storeOrig XC --poolSize 0 --noph 2> results/testClipOverlapCoordPool0ClipStatsClipsOnly.log
if [ $? != 8 ]
then
    status=1
    echo did not get expected return value for no pool.
fi
diff results/testClipOverlapCoordPool0ClipStatsClipsOnly.sam expected/testClipOverlapCoordPool0Clip.sam
let "status |= $?"
diff results/testClipOverlapCoordPool0ClipStatsClipsOnly.log expected/testClipOverlapCoordPool0ClipStats.log
let "status |= $?"

# Test clipping files sorted by read name.
../bin/bam clipOverlap --stats --readName --in testFiles/testClipOverlapReadName1.sam --out results/testClipOverlapReadNameStats1.sam --storeOrig XC --noph 2> results/testClipOverlapReadNameStats1.log
let "status |= $?"
diff results/testClipOverlapReadNameStats1.sam expected/testClipOverlapReadName1.sam
let "status |= $?"
diff results/testClipOverlapReadNameStats1.log expected/testClipOverlapReadNameStats.log
let "status |= $?"

# Test clipping files sorted by read name.
../bin/bam clipOverlap --stats --readName --in testFiles/testClipOverlapReadName1.sam --out results/testClipOverlapReadNameStats1Ex0.sam --storeOrig XC --excludeFlag 0 --noph 2> results/testClipOverlapReadNameStats1Ex0.log
let "status |= $?"
diff results/testClipOverlapReadNameStats1Ex0.sam expected/testClipOverlapReadName1Ex0.sam
let "status |= $?"
diff results/testClipOverlapReadNameStats1Ex0.log expected/testClipOverlapReadNameStatsEx0.log
let "status |= $?"

# Test clipping files sorted by coordinate
../bin/bam clipOverlap --stats --in testFiles/testClipOverlapCoord1.sam --out results/testClipOverlapCoordStats1.sam --storeOrig XC --noph 2> results/testClipOverlapCoordStats1.log
let "status |= $?"
diff results/testClipOverlapCoordStats1.sam expected/testClipOverlapCoord1.sam
let "status |= $?"
diff results/testClipOverlapCoordStats1.log expected/testClipOverlapCoordStats.log
let "status |= $?"

../bin/bam clipOverlap --stats --in testFiles/testClipOverlapCoord1.sam --out results/testClipOverlapCoordStats1Ex0.sam --storeOrig XC --excludeFlag 0x0000 --noph 2> results/testClipOverlapCoordStats1Ex0.log
let "status |= $?"
diff results/testClipOverlapCoordStats1Ex0.sam expected/testClipOverlapCoord1Ex0.sam
let "status |= $?"
diff results/testClipOverlapCoordStats1Ex0.log expected/testClipOverlapCoordStatsEx0.log
let "status |= $?"

if [ $status != 0 ]
then
  echo failed testClipOverlap.sh
  exit 1
fi

