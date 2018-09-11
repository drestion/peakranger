/*
 * bamFormatAux.cpp
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#include "bamFormatAux.h"

using namespace std;
using namespace BamTools;

namespace reads {

string getR1Chr(const BamAlignment& read, const RefVector& refvec) {
    return getRChr(read.RefID, refvec);
}

string getRChr(const int32_t& ref_id, const RefVector& refvec) {
    if (ref_id < (int32_t) refvec.size() && ref_id >= 0) {
        return refvec[ref_id].RefName;
    }
    return string("");
}

bool isProperPERead(const BamAlignment& read) {
    return isFirstPEGoodRead(read) && read.IsProperPair();
}

bool isGoodSERead(const BamAlignment& read) {
    return isSERead(read) && !(read.IsFailedQC()) && !(read.IsDuplicate())
            && read.IsPrimaryAlignment();
}

bool bothArePosReads(const bool& dir, const bool& dir2) {
    return dir && dir2;
}

bool bothAreNegReads(bool& dir, bool& dir2) {
    return !dir && !dir2;
}

void processAbnormalRead(const BamAlignment& read, const RefVector& ref) {
    return;
}

Strand getR1Strand(const BamTools::BamAlignment& read) {
    return Strand(!(read.IsReverseStrand()));
}

Strand getR2Strand(const BamTools::BamAlignment& read) {
    return Strand(!(read.IsMateReverseStrand()));
}

bool isChrToParse(vector<string>& chrs_to_parse, string& chr) {
    return std::find(chrs_to_parse.begin(), chrs_to_parse.end(), chr)
            != chrs_to_parse.end();
}

bool isFirstPEGoodRead(const BamAlignment& read) {
    return isGoodPERead(read) && read.IsFirstMate();
}

string getR2Chr(const BamAlignment& read, const RefVector& refvec) {
    return getRChr(read.MateRefID, refvec);
}

void getBamHitsFromPEBamRead(const BamAlignment& read, const RefVector& ref,
        BamHit& bam1, BamHit& bam2) {
    bam1 = getRead1FromPEBamRead(read, ref);
    bam2 = getRead2FromPEBamRead(read, ref);
}

void getR1StartEnd(const BamTools::BamAlignment& read, int32_t& start,
        int32_t& end) {
    start = getR1Start(read);
    end = getR1End(read);
}

int32_t getR1Start(const BamTools::BamAlignment& read) {
    return read.Position;
}

int32_t getR1End(const BamTools::BamAlignment& read) {
    return read.GetEndPosition();
}

BamHit getRead1FromPEBamRead(const BamTools::BamAlignment& read,
        const BamTools::RefVector& ref) {
    return BamHit(read, ref);
}

int32_t getR2Start(const BamTools::BamAlignment& read) {
    return read.MatePosition;
}

int32_t getR2End(const BamTools::BamAlignment& read) {
    return guessMateEnd(read);
}

void getR2StartEnd(const BamTools::BamAlignment& read, int32_t& start,
        int32_t& end) {
    start = getR2Start(read);
    end = getR2End(read);
}

int32_t guessMateEnd(const BamTools::BamAlignment& read) {
    return read.MatePosition + getReadLength(read);
}

BamHit getRead2FromPEBamRead(const BamTools::BamAlignment& read,
        const BamTools::RefVector& ref) {
    BamHit hit;
    hit.setChr(getR2Chr(read, ref));
    hit.setStart(getR2Start(read));
    hit.setEnd(getR2End(read));
    hit.setStrand(getR2Strand(read));
    return hit;
}

bool isSecondPEGoodRead(const BamTools::BamAlignment& read) {
    return isGoodPERead(read) && !read.IsFirstMate();
}

bool isGoodPERead(const BamTools::BamAlignment& read) {
    return isPERead(read) && read.IsMapped() && !(read.IsFailedQC())
            && !(read.IsDuplicate()) && (read.IsPrimaryAlignment());
}

bool isSameReadName(const BamTools::BamAlignment& lhs,
        const BamTools::BamAlignment& rhs) {
    return lhs.Name == rhs.Name;
}

bool isPERead(const BamTools::BamAlignment& read) {
    return read.IsPaired();
}

bool isSERead(const BamTools::BamAlignment& read) {
    return !isPERead(read);
}

bool isGoodRead(const BamTools::BamAlignment& read) {
    return read.IsMapped() && !(read.IsFailedQC()) && !(read.IsDuplicate())
            && read.IsPrimaryAlignment();
}

bool isDifferentReadName(const BamTools::BamAlignment& lhs,
        const BamTools::BamAlignment& rhs) {
    return !isSameReadName(lhs, rhs);
}

}

uint32_t reads::getReadLength(const BamTools::BamAlignment& read) {
	//not end - pos + 1; because bamtools's end is end-1
	//HWUSI-EAS1758_6:1:12:4982:11510:0:1:1	0	chr11	189734	255	36M	*	0	0	*	*
	//will be 189733 instead of 189734
	return read.GetEndPosition() - read.Position;

}
