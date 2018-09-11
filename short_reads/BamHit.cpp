/*
 * BamHit.cpp
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#include "BamHit.h"
#include "common/stl_header.h"
#include "bamFormatAux.h"
using namespace std;

namespace reads {

BamHit::BamHit() {

}

BamHit::BamHit(const BamTools::BamAlignment& bam, const BamTools::RefVector& ref) {
    getStartEndDir(bam);
    mChr = getR1Chr(bam, ref);
}

BamHit::~BamHit() {

}

BamHit::BamHit(const BamTools::BamAlignment& bam, const std::string& chr) {
    getStartEndDir(bam);
    mChr = chr;
}

BamHit::BamHit(const BamTools::BamAlignment& bam, const char* chr) {
    getStartEndDir(bam);
    mChr = string(chr);
}

int32_t BamHit::getReadLength(const BamTools::BamAlignment& bam) {
    // Using Cigar string to get the end position
    return reads::getReadLength(bam);
}

void BamHit::getStartEndDir(const BamTools::BamAlignment& bam) {
    mStart = bam.Position;
    mEnd = getReadLength(bam) + mStart;
    mStrand = !(bam.IsReverseStrand());
}
/* namespace bams */

}

