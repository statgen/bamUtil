../bin/bam clipOverlap --in testFiles/testClipOverlap.sam --out results/testClipOverlap.sam --storeOrig XC 2> results/testClipOverlap.log && diff results/testClipOverlap.sam expected/testClipOverlap.sam && diff results/testClipOverlap.log expected/testClipOverlap.log


