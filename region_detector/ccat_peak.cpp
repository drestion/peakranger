/*
 * ccat_peakt.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#include "ccat_peak.h"
#include <tab_file/Region.h>

namespace ccat_aux {

peak_t::peak_t() :
        chromIndex(0), start(0), end(0), peak(0), l1Count(0), l2Count(0), reSampledL1Count(0), reSampledL2Count(0), foldChange(
                0), qValue(0), isSignificant(0) {

}

peak_t::~peak_t() {

}

bool lessFoldChange(const peak_t& lhs, const peak_t& rhs) {
    return lhs.foldChange < rhs.foldChange;
}

std::ostream& operator<<(std::ostream& os, peak_t& pk) {
    tab_file::Region region(pk.start, pk.end);
    os << region << "\t" << pk.foldChange << "\t" << pk.qValue << "\t" << pk.isSignificant;
    return os;
}

bool greaterFoldChange(const peak_t & lhs, const peak_t & rhs) {
    return lhs.foldChange > rhs.foldChange;
}

} /* namespace ccat_aux */
