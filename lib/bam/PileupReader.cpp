#include <stdexcept>
/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "PileupReader.h"

bool ReadBase::set(SamRecord &samRecord, int basePosition)
{
    throw std::logic_error("ReadBase::set unimplemented");
    _samRecord = &samRecord;
    return false;
}

bool SequenceCoverageReader::getBases(
    int referenceID,
    uint32_t position,
    std::vector<ReadBase> &readBases,
    bool &basesInserted)
{
    throw std::logic_error("SequenceCoverageReader::getBases unimplemented");
    return true;
}

bool SequenceCoverageReader::getInsertedBases(
    int referenceID,
    uint32_t position,
    std::vector<ReadInsertion> &readInsertions
    )
{
    abort();
    return true;
}
