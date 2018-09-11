/*
 * Read.h
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#ifndef READ_H_
#define READ_H_

#include "common/stl_header.h"
#include "short_reads/Strand.h"
namespace reads {

class Read {
	friend std::ostream& operator<<(std::ostream& os, const Read& rhs) {

		os << rhs.mChr << "\t" << rhs.mStart << "\t" << rhs.mEnd << "\t"
				<< rhs.mStrand << "\n";
		return os;
	}
	friend bool operator!=(const Read& lhs, const Read& rhs) {
		return !lhs.operator ==(rhs);
	}

public:
	Read();
	Read(const int32_t& start, const int32_t& end, const char* chr,
			const Strand& strand);
	virtual ~Read();

	void set(const int32_t& start, const int32_t& end, const char* chr,
			const Strand& strand);

	Strand getStrand() const;

	void setStrand(Strand mStrand);

	int32_t getEnd() const {
		return mEnd;
	}

	int32_t getStart() const {
		return mStart;
	}
	std::string getChr() const {
		return mChr;
	}

	int32_t getLength() const {
		return mEnd - mStart;
	}
	void setChr(std::string chr) {
		mChr = chr;
	}

	void setEnd(int32_t end) {
		mEnd = end;
	}

	void setStart(int32_t start) {
		mStart = start;
	}
	bool operator==(const Read& rhs) const {
		return mStart == rhs.mStart &&
		       mEnd == rhs.mEnd &&
		       mChr == rhs.mChr &&
		       mStrand == rhs.mStrand;
	}

protected:
	int32_t mStart;
	int32_t mEnd;
	Strand mStrand;
	std::string mChr;
};

}

/* namespace reads */
#endif /* READ_H_ */
