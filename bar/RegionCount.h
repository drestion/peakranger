/*
 * RegionCount.h
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#ifndef REGIONCOUNT_H_
#define REGIONCOUNT_H_
#include "concepts/RegionInt32.h"

#include <ostream>
namespace ranger {
namespace bar {

class RegionCount {
	friend std::ostream& operator<<(std::ostream& os,  RegionCount& cnt) {
		os << cnt.mRegion << "\t" << cnt.mCnt;
		return os;
	}
public:
	RegionCount();
	RegionCount(const concepts::RegionInt32& reg, const int32_t& cnt);
	virtual ~RegionCount();

	int32_t getCnt() const {
		return mCnt;
	}

	void setCnt(int32_t cnt) {
		mCnt = cnt;
	}

	concepts::RegionInt32 getRegion() const {
		return mRegion;
	}

	void setRegion(const concepts::RegionInt32& region) {
		mRegion = region;
	}

private:
	concepts::RegionInt32 mRegion;
	int32_t mCnt;
};

} /* namespace bar */
} /* namespace ranger */
#endif /* REGIONCOUNT_H_ */
