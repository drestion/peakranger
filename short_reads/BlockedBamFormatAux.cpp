/*
 * BlockedBamFormatAux.cpp
 *
 *  Created on: May 11, 2012
 *      Author: xin
 */

#include "BlockedBamFormatAux.h"
#include "bamFormatAux.h"
#include "common/ranger_debug.h"
#include <algorithm>
using namespace reads;
using namespace BamTools;
using namespace std;
using namespace boost;
using namespace ranger::concepts;

string reads::getBlockedR1Chr(const BamAlignment & read,
		const RefVector & refvec) {
}

string reads::getBlockedR2Chr(const BamAlignment & read,
		const RefVector & refvec) {
}

string reads::getBlockedRChr(const int32_t & ref_id, const RefVector & refvec) {
}

int32_t reads::getBlockedR1Start(const BamAlignment & read) {
}

int32_t reads::getBlockedR1End(const BamAlignment & read) {
}

void reads::getBlockedR1StartEnd(const BamAlignment & read, int32_t & start,
		int32_t & end) {
}

int32_t reads::getBlockedR2Start(const BamAlignment & read) {
}

int32_t reads::getBlockedR2End(const BamAlignment & read) {
}

void reads::getBlockedR2StartEnd(const BamAlignment & read, int32_t & start,
		int32_t & end) {
}

Strand reads::getBlockedR1Strand(const BamAlignment & read) {
}

Strand reads::getBlockedR2Strand(const BamAlignment & read) {
}

BlockedBamHit reads::getBlockedRead1FromPEBamRead(const BamAlignment & read,
		const RefVector & ref) {
	return BlockedBamHit(read, ref);
}

BlockedBamHit reads::getBlockedRead2FromPEBamRead(const BamAlignment & read,
		const RefVector & ref) {
  cerr << "Not implemented!\n";
  
  return BlockedBamHit(read, ref);
}

void reads::getBlockedBamHitsFromPEBamRead(const BamAlignment & read,
		const RefVector & ref, BlockedBamHit & r1, BlockedBamHit & r2) {
	r1 = getBlockedRead1FromPEBamRead(read, ref);
	r2 = getBlockedRead2FromPEBamRead(read, ref);
}

RegionInt32 reads::regionFromCigar(const CigarString& cigar) {
	return cigar.toRegion();
}

Read reads::readFromCigar(const CigarString& cigar) {

	return readFromCigar(cigar, "", Strand(true));
}

void reads::getCigarString(const BamAlignment& bam,
		vector<boost::shared_ptr<CigarString> >& lica) {
	int32_t newOffset = bam.Position;
	foreach(CigarOp c, bam.CigarData) {
		lica.push_back(CigarStringFactory::buildCigar(c, newOffset));
		if (c.Type == 'M' || c.Type == 'D' || c.Type == 'N') {
			newOffset += c.Length;
		}
	}
}

bool reads::isAdditiveCigar(const CigarString& c) {
//	cout <<"[isAdditiveCigar]in it\n";
	RegionInt32 r = c.toRegion();
	return r.getL() != r.getR();
}

Read reads::readFromCigar(const CigarString& cigar, const string& chr,
		const Strand& strand) {
	RegionInt32 r = regionFromCigar(cigar);
//	cout <<"[readFromCigar]after regionFromCigar\n";
	Read rd(r.getL(), r.getR(), chr.c_str(), strand);
	return rd;
}

void reads::parseBamBlocks(vector<Read>& resBlocks, const BamAlignment& bam,
		const RefVector & ref) {
	vector<boost::shared_ptr<CigarString> > lica;

	getCigarString(bam, lica);
//	cout <<"[parseBamBlocks]after getCigarString\n";
	for (size_t i = 0; i < lica.size(); i++) {

		if (isAdditiveCigar(*lica[i])) {
//			cout <<"[parseBamBlocks]after getCigarString\n";
			resBlocks.push_back(
					readFromCigar(*lica[i], getR1Chr(bam, ref),
							!bam.IsReverseStrand()));
		}
	}
}

void reads::lessGoFirst(BamAlignment& mread, BamAlignment& read) {
	if (mread.Position < read.Position) {
		std::swap(read, mread);
	}
}

std::ostream& BamTools::operator <<(std::ostream& os, const BamAlignment& bam) {
	if (bam.IsMapped()) {
		os << bam.Position << "\t" << bam.Name;
	} else {
		os << bam.Name;
	}
	return os;
}

void reads::extendRead(BamTools::BamAlignment& r, int32_t length) {
	if (length > 0) {
		r.Length += (length - r.Length);
	} else {
		r.Length += length;
	}
	r.Length = r.Length > 0 ? r.Length : 0;
//    if (r.CigarData.size()) {
//        CigarOp& lastOp(r.CigarData.back());
//        lastOp.Length += (length - lastOp.Length);
//    }
//    BamTools::RefVector foo;
//    BlockedBamHit hit(r, foo);
//    vector<Read> rds;
//    hit.getBlocks(rds);
//    if (rds.size()) {
//        int32_t start = rds.front().getStart();
//        int32_t end = rds.back().getEnd();
//        if (length > 0) {
//            end += (length - (end - start));
//        }
//        end = end > 0 ? end : 0;
//        start = end > 0 ? start : 0;
//        rds.back().setEnd(end);
//        rds.front().setStart(start);
//
//    }
}

void reads::extendBamRead(const BamTools::BamAlignment& read, vector<Read>& rds,
		uint32_t length) {
	BamTools::RefVector foo;
	rds.resize(0);
	parseBamBlocks(rds,read,foo);

	if (rds.size()) {
		int32_t start = rds.front().getStart();
		int32_t end = rds.back().getEnd();
		//TODO: the logic here is not the same as in
		// FDR_THRESHOLDER's _cal_p_read_end
		if (read.IsReverseStrand()) {
			if (end - start < (int32_t)length) {
				start -= (length - (end - start));
			}
			start = start > 0 ? start : 0;
			rds.front().setStart(start);
		} else {
			if (end - start < (int32_t)length) {
				end += (length - (end - start));
			}
			rds.back().setEnd(end);
		}
	}
}
