#ifndef _MAPPERUSEROPTION_H_
#define _MAPPERUSEROPTION_H_

struct MapperUserOption
{
    MapperUserOption() :
            checkIndexWordMutations(false),
            debug(false),
            expectedSNPRate(.001),
            expectedErrorRate(.01),
            expectedInDelRate(.0001),
            mismatchCutoff(-1),
            minimumIndexCount(2),
            maximumIndexCount(4),
            forceSmithWaterman(false),
            allowSmithWaterman(false),
            smithWatermanBandSize(4),
            readSumQCutoff(-1),
            showReferenceBases(false),
            qValueCutoff(INT32_MAX),
            genomePositionFilterWidth(1500),
            trimLeft(0),
            trimRight(0),
            qualityTrim(0)
    {
    }
    bool     checkIndexWordMutations;
    bool     debug;

    double   expectedSNPRate;
    double   expectedErrorRate;
    double   expectedInDelRate;
    // double   expectedInDelErrorRate; // it is never used

    int      mismatchCutoff;
    uint32_t minimumIndexCount;
    uint32_t maximumIndexCount;

    bool     forceSmithWaterman;
    bool     allowSmithWaterman;
    int      smithWatermanBandSize;
    int         readSumQCutoff;
    bool     showReferenceBases;
    int32_t  qValueCutoff;

    uint32_t   genomePositionFilterWidth;

    int trimLeft, trimRight;

    int qualityTrim;

    std::string readGroupID;
};


#endif /* _MAPPERUSEROPTION_H_ */
