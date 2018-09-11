/*
 * RegionCount.cpp
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#include "bar/RegionCount.h"

namespace ranger {
namespace bar {

RegionCount::RegionCount() {
}

RegionCount::RegionCount(const concepts::RegionInt32& reg, const int32_t& cnt) :
		mRegion(reg), mCnt(cnt) {
}

RegionCount::~RegionCount() {
}

} /* namespace bar */
} /* namespace ranger */
