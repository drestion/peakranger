/*
 * BamBlockPEReadsParser_hotelImphotelImp.cpp
 *
 *  Created on: May 13, 2012
 *      Author: xfeng
 */

#include "BamBlockPEReadsParserhotelImp.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "short_reads/BlockedBamFormatAux.h"
#include "short_reads/bamFormatAux.h"
#include "short_reads/BamHit.h"
#include "short_reads/ReadPair.h"
#include "short_reads/BlockedBamHit.h"
#include "BamParserAux.h"
#include "common/ranger_debug.h"
#include <algorithm>

using namespace reads;
using namespace BamTools;
using namespace std;
namespace parser {
namespace aux {

BamBlockPEReadsParser_hotelImp::BamBlockPEReadsParser_hotelImp() {

}

BamBlockPEReadsParser_hotelImp::~BamBlockPEReadsParser_hotelImp() {
}

bool BamBlockPEReadsParser_hotelImp::checkedIn(const BamAlignment & read) {
	void adjustPosition(BamAlignment& mread, const BamAlignment& read);
	foreach(BamAlignment & b, mUnMatched) {
		if (isSameReadName(read, b)) {
			return true;
		}
	}
	return false;
}

void BamBlockPEReadsParser_hotelImp::checkIn(
		const BamTools::BamAlignment & bam) {
//	cout << "Check in:" << (int)bam.Position << "\t" << bam.Name << "\n";
	mUnMatched.push_back(bam);
}

void BamBlockPEReadsParser_hotelImp::parse(const BamAlignment & read,
		PairEndedReads<BlockedRead> & reads, const RefVector & ref) {
	if (isGoodPERead(read)) {
		if (!checkedIn(read)) {
			checkIn(read);
		} else {
			BamAlignment mread = checkOut(read);
			if (mread.Position < read.Position) {
				insertRead(mread, read, reads, ref);
			} else {
				insertRead(read, mread, reads, ref);
			}
		}
	} else if (isGoodSERead(read)) {
		throw NotPairEndBamRead("The bam file contains single end reads.");
	} else {
		processAbnormalRead(read, ref);
	}
}

void BamBlockPEReadsParser_hotelImp::flush(
		reads::PairEndedReads<reads::BlockedRead>& reads,
		const BamTools::RefVector& ref) {
	foreach(BamAlignment& b, mUnMatched) {
//		cout << "flushed :" << (int)b.Position << "\t" << b.Name << "\n";
		insertRead(b, BamAlignment(), reads, ref);
	}
}

BamTools::BamAlignment BamBlockPEReadsParser_hotelImp::checkOut(
		const BamTools::BamAlignment & bam) {
	BamAlignment res;
	vector<BamAlignment>::iterator it;
	it = mUnMatched.begin();
	while (it != mUnMatched.end()) {
		if (isSameReadName(bam, *it)) {
			res = *it;
			break;
		}
		++it;
	}

	if (it != mUnMatched.end()) {
//		cout << "found the mate for:" << (int)bam.Position << "\t" << bam.Name
//				<< "\tas\t" << it->Position<<"\t"<<it->Name<<"\n";
		mUnMatched.erase(it);
	} else {
//		cout << "didnt find the mate for:" << (int)bam.Position << "\t" << bam.Name
//				<< "\n";
	}
	return res;
}
} /* namespace aux */
} /* namespace parser */
