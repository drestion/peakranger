/*
 * Read.cpp
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#include "Read.h"

using namespace reads;
using namespace std;
reads::Read::Read() :
		mStart(0), mEnd(0), mStrand(true), mChr("") {
}

reads::Read::~Read() {
}

Strand reads::Read::getStrand() const {
	return mStrand;
}

reads::Read::Read(const int32_t & start, const int32_t & end, const char *chr,
		const Strand & strand) :
		mStart(start), mEnd(end), mStrand(strand), mChr(string(chr)) {
}

void reads::Read::set(const int32_t& start, const int32_t& end, const char* chr,
		const Strand& strand) {
	setStart(start);
	setEnd(end);
	setChr(string(chr));
	setStrand(strand);
}

void reads::Read::setStrand(Strand mStrand) {
	this->mStrand = mStrand;
}

