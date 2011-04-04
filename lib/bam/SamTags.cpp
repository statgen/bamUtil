/*
 *  Copyright (C) 2011  Regents of the University of Michigan
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

#include "SamTags.h"

const char* SamTags::BQ_TAG = "BQ";
const char SamTags::BQ_TAG_TYPE = 'Z';
const char* SamTags::ORIG_POS_TAG = "OP";
const char SamTags::ORIG_POS_TAG_TYPE = 'i';
const char* SamTags::ORIG_CIGAR_TAG = "OC";
const char SamTags::ORIG_CIGAR_TAG_TYPE = 'Z';
const char* SamTags::ORIG_QUAL_TAG = "OQ";
const char SamTags::ORIG_QUAL_TAG_TYPE = 'Z';
