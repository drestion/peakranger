/*
 * ccat_profile_helper.h
 *
 *  Created on: Apr 5, 2012
 *      Author: xin
 */

#ifndef CCAT_PROFILE_HELPER_H_
#define CCAT_PROFILE_HELPER_H_

#include "ccat_profile.h"
#include "common/stl_header.h"
#include "chrt.h"
namespace ccat_aux {

void buildL1Profile(ccat_profile_t& profile, const chr_t& chr, std::vector<size_t>& result);

void buildL2Profile(ccat_profile_t& profile, const chr_t& chr, std::vector<size_t>& result);

void buildrsL2Profile(ccat_profile_t& profile, const chr_t& chr, const double ratio, std::vector<size_t>& result);

void buildrsL1Profile(ccat_profile_t& profile, const chr_t& chr, const double ratio, std::vector<size_t>& result);

bool rsRatio(const double ratio);

template<typename AT>
void getMaxProfilePnt(size_t& max, AT& vec) {
    size_t t1;
    t1 = *std::max_element(vec.begin(), vec.end());

    max = t1 > max ? t1 : max;
}

}

#endif /* CCAT_PROFILE_HELPER_H_ */
