#include <cassert>
#include "Chromosome.h"

Chromosome::Chromosome(GenomeSequence* gs, uint chromosomeIndex)
{
    assert(gs);
    assert(chromosomeIndex < (uint)gs->getChromosomeCount());

    this->gs = gs;
    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}

Chromosome::Chromosome(GenomeSequence* gs, const char* chromosomeName)
{
    assert(gs);
    this->gs = gs;

    int chromosomeIndex = gs->getChromosome(chromosomeName);
    assert(chromosomeIndex != INVALID_CHROMOSOME_INDEX);

    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}

Chromosome::Chromosome(const char* genomseSequenceFileName, uint chromosomeIndex, bool isColorSpace) 
{
    std::string s(genomseSequenceFileName);
    if (this->gs) delete gs;
    gs = new GenomeSequence;
    assert(gs);
    gs->setReferenceName(s);
    assert(!gs->open(isColorSpace));

    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}

Chromosome::Chromosome(const std::string& genomseSequenceFileName, uint chromosomeIndex, bool isColorSpace) 
{
    if (this->gs) delete gs;
    gs = new GenomeSequence;
    assert(gs);
    gs->setReferenceName(genomseSequenceFileName);
    assert(!gs->open(isColorSpace));

    this->offset = gs->getChromosomeStart((int)chromosomeIndex);
    this->chromosomeSize = gs->getChromosomeSize((int)chromosomeIndex);
}
