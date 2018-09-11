/*
 * CCATPeakFinder.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xfeng
 */

#ifndef CCATPEAKFINDER_H_
#define CCATPEAKFINDER_H_

#include "common/stl_header.h"
#include "chrt.h"
#include "ccat_config.h"
#include "ccat_peak.h"
#include "ccat_profile.h"

namespace ccat_aux {

class CCATPeakFinder {
public:
    CCATPeakFinder();
    virtual ~CCATPeakFinder();

    //PeakFinding: peak processing, will be called by the main routine of CCAT
    int PeakFinding(std::vector<ccat_aux::chr_t>&chroms, size_t chromNum, double l1Ratio, double l2Ratio,
            size_t &maxL1Count, size_t &maxL2Count, const ccat_aux::ccat_config_t& config);

private:
    //GetPeaksInOneChrom1: strand-insensitive mode: get peak location from tags for one chromosome
    int GetPeaksInOneChrom1(ccat_aux::chr_t& chrom, double l1Ratio, double l2Ratio, size_t &maxL1Count,
            size_t &maxL2Count, const ccat_aux::ccat_config_t& config);

    //GetLocalMaxima: find the local maxima in the profile
    int GetLocalMaxima(const std::vector<size_t>& profile, std::vector<ccat_aux::peak_t>& peaks, const int minDist,
            const size_t minCount);

    void callPeaks(const std::vector<size_t>& rsProfile1, const std::vector<size_t>& rsProfile2,
            const std::vector<size_t>& profile1, const std::vector<size_t>& profile2, const ccat_config_t& config,
            const size_t chromSize, std::vector<ccat_aux::peak_t>& result);

    bool hasLargerNeighbors(size_t minDist, int& tmpStart, int& tmpEnd, const std::vector<size_t>& profile,
            const peak_t& pk);

    ccat_profile_t profile;

};

} /* namespace ccat_aux */
#endif /* CCATPEAKFINDER_H_ */
