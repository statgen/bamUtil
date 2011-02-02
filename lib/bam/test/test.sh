./samTest 2> results/samTest.log && \
diff expected/testEqWithBases.sam results/outSamEqBases.sam && \
diff expected/testEqWithEq.sam results/outSamEqEquals.sam && \
diff expected/testEqWithOrig.sam results/outSamEqOrig.sam && \
diff expected/testEqWithBases.bam results/outSamEqBases.bam && \
diff expected/testEqWithEq.bam results/outSamEqEquals.bam && \
diff expected/testEqWithOrig.bam results/outSamEqOrig.bam && \
diff expected/testEqWithBases.sam results/outBamEqBases.sam && \
diff expected/testEqWithEq.sam results/outBamEqEquals.sam && \
diff expected/testEqWithOrig.sam results/outBamEqOrig.sam && \
diff expected/testEqWithBases.bam results/outBamEqBases.bam && \
diff expected/testEqWithEq.bam results/outBamEqEquals.bam && \
diff expected/testEqWithOrig.bam results/outBamEqOrig.bam && \
diff expected/samTest.log results/samTest.log
