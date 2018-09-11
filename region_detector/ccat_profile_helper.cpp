/*
 * ccat_profile_helper.cpp
 *
 *  Created on: Apr 5, 2012
 *      Author: xin
 */

#include "ccat_profile_helper.h"

namespace ccat_aux {

bool rsRatio(const double ratio) {
    return rand() > RAND_MAX * ratio;
}

void buildL1Profile(ccat_profile_t & profile, const chr_t & chr, std::vector<size_t>& result) {

    profile.GetProfile(chr.l1PosTags, result, result.size());
    profile.GetProfile(chr.l1NegTags, result, result.size(), true);
}

void buildL2Profile(ccat_profile_t & profile, const chr_t & chr, std::vector<size_t>& result) {

    profile.GetProfile(chr.l2PosTags, result, result.size(), false);
    profile.GetProfile(chr.l2NegTags, result, result.size(), true);
}

void buildrsL2Profile(ccat_profile_t & profile, const chr_t & chr, const double ratio, std::vector<size_t> & result) {
    profile.ResampleProfile(chr.l2PosTags, result, result.size(), ratio, false);
    profile.ResampleProfile(chr.l2NegTags, result, result.size(), ratio, true);
}

void buildrsL1Profile(ccat_profile_t & profile, const chr_t & chr, const double ratio, std::vector<size_t> & result) {

    profile.ResampleProfile(chr.l1PosTags, result, result.size(), ratio, false);
    profile.ResampleProfile(chr.l1NegTags, result, result.size(), ratio, true);
}

}
