cmake_minimum_required(VERSION 3.5)
project(bamUtil)

include_directories(${CMAKE_SOURCE_DIR}/../libStatGen/include)

execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)
set(VERSION "1.0.14")

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DVERSION='\"${VERSION}\"' -DUSER='\"${USER}\"' -DDATE='\"${DATE}\"' -D__ZLIB_AVAILABLE__")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(SOURCE_FILES
        src/Bam2FastQ.cpp
        src/Bam2FastQ.h
        src/BamExecutable.cpp
        src/BamExecutable.h
        src/ClipOverlap.cpp
        src/ClipOverlap.h
        src/Convert.cpp
        src/Convert.h
        src/Covariates.h
        src/Dedup.cpp
        src/Dedup.h
        src/Dedup_LowMem.cpp
        src/Dedup_LowMem.h
        src/Diff.cpp
        src/Diff.h
        src/DumpHeader.cpp
        src/DumpHeader.h
        src/DumpIndex.cpp
        src/DumpIndex.h
        src/DumpRefInfo.cpp
        src/DumpRefInfo.h
        src/ExplainFlags.cpp
        src/ExplainFlags.h
        src/Filter.cpp
        src/Filter.h
        src/FindCigars.cpp
        src/FindCigars.h
        src/GapInfo.cpp
        src/GapInfo.h
        src/HashErrorModel.cpp
        src/HashErrorModel.h
        src/Logger.cpp
        src/Logger.h
        src/LogisticRegression.cpp
        src/LogisticRegression.h
        src/Main.cpp
        src/MateMapByCoord.cpp
        src/MateMapByCoord.h
        src/MathCholesky.cpp
        src/MathCholesky.h
        src/MergeBam.cpp
        src/MergeBam.h
        src/OverlapClipLowerBaseQual.cpp
        src/OverlapClipLowerBaseQual.h
        src/OverlapHandler.cpp
        src/OverlapHandler.h
        src/PileupElementBaseQCStats.cpp
        src/PileupElementBaseQCStats.h
        src/PolishBam.cpp
        src/PolishBam.h
        src/Prediction.cpp
        src/Prediction.h
        src/ReadIndexedBam.cpp
        src/ReadIndexedBam.h
        src/ReadReference.cpp
        src/ReadReference.h
        src/Recab.cpp
        src/Recab.h
        src/Revert.cpp
        src/Revert.h
        src/SplitBam.cpp
        src/SplitBam.h
        src/SplitChromosome.cpp
        src/SplitChromosome.h
        src/Squeeze.cpp
        src/Squeeze.h
        src/Stats.cpp
        src/Stats.h
        src/TrimBam.cpp
        src/TrimBam.h
        src/Validate.cpp
        src/Validate.h
        src/WriteRegion.cpp
        src/WriteRegion.h)

add_executable(bam ${SOURCE_FILES})
target_link_libraries(bam ${CMAKE_SOURCE_DIR}/../libStatGen/libStatGen.a)
target_link_libraries(bam libz.dylib)