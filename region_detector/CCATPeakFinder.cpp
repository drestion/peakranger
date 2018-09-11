/*
 * CCATPeakFinder.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xfeng
 */

#include "CCATPeakFinder.h"
#include "common/ranger_debug.h"
#include "ccat_profile_helper.h"
using namespace std;
namespace ccat_aux {

CCATPeakFinder::CCATPeakFinder() {

}

CCATPeakFinder::~CCATPeakFinder() {

}

int CCATPeakFinder::PeakFinding(vector<chr_t>& chroms, size_t chromNum,
        double l1Ratio, double l2Ratio, size_t& maxL1Count, size_t& maxL2Count,
        const ccat_config_t& config) {
    MARK_FUN("CCATPeakFinder::PeakFinding");

    size_t i;
    size_t peakNum, tmpI;

    maxL1Count = 0;
    maxL2Count = 0;
    peakNum = 0;

    for (i = 0; i < chromNum; i++) {

        if (config.isStrandSensitiveMode) {

            cout << "TF mode hasnt been implemented\n";
            throw string("TF bad");
        } else {

            tmpI = GetPeaksInOneChrom1(chroms[i], l1Ratio, l2Ratio, maxL1Count,
                    maxL2Count, config);
        }

        //TODO: Check if this is necessary
        if (tmpI < 0) {
            return -1;
        }

        cout << "Discovered " << tmpI << " candidate peaks in "
                << chroms[i].chromName << "\n";

        peakNum += tmpI;
    }

    return peakNum;
}

void CCATPeakFinder::callPeaks(const vector<size_t>& rsProfile1,
        const vector<size_t>& rsProfile2, const vector<size_t>& profile1,
        const vector<size_t>& profile2, const ccat_config_t& config,
        const size_t chromSize, vector<peak_t>& result) {
    MARK_FUN("CCATPeakFinder::callPeaks");
    result.resize(0);
    size_t minDist = config.slidingWinSize / config.movingStep + 1;
    size_t maxPossiblePeaks = chromSize / config.movingStep + 1;

    vector<peak_t> tmp(maxPossiblePeaks);
    size_t tmpPeakNum = GetLocalMaxima(rsProfile1, tmp, minDist,
            config.minCount);

    assert(tmpPeakNum <= maxPossiblePeaks);
    LOG_DEBUG1("tmpPeakNum "<<tmpPeakNum);
    foreach(peak_t pk, tmp) {
        pk.l1Count = profile1[pk.peak];
        pk.l2Count = profile2[pk.peak];
        pk.reSampledL1Count = rsProfile1[pk.peak];
        pk.reSampledL2Count = rsProfile2[pk.peak];
        if ((pk.reSampledL1Count >= config.minCount)) {
            result.push_back(pk);
        }
    }
}

//GetPeaksInOneChrom1: strand-insensitive mode: get peak location from tags for one chromosome
int CCATPeakFinder::GetPeaksInOneChrom1(chr_t& chrom, double l1Ratio,
        double l2Ratio, size_t& maxL1Count, size_t& maxL2Count,
        const ccat_config_t& config) {

    MARK_FUN("Entering ccat::GetPeaksInOneChrom1");
    size_t profileLen = (chrom.chromSize) / config.movingStep + 1;
    LOG_DEBUG1("profileLen "<<profileLen);LOG_DEBUG1("config.fragmentSize "<<config.fragmentSize);
    vector<size_t> profile1(profileLen, 0);
    vector<size_t> profile2(profileLen, 0);
    vector<size_t> rsProfile1(profileLen, 0);
    vector<size_t> rsProfile2(profileLen, 0);

    size_t maxPossiblePeaks = (chrom.chromSize) / config.movingStep + 1;
    vector<peak_t> tmpPeaks(maxPossiblePeaks);
    size_t minDist;

    profile.setFragmentSize(config.fragmentSize);
    profile.setMovingStep(config.movingStep);
    profile.setProfileLength(profileLen);

    buildL1Profile(profile, chrom, profile1);
    buildL2Profile(profile, chrom, profile2);
    buildrsL1Profile(profile, chrom, l1Ratio, rsProfile1);
    buildrsL2Profile(profile, chrom, l2Ratio, rsProfile2);

    profile.SmoothProfile(profile1, profileLen,
            config.slidingWinSize / config.movingStep / 2);
    profile.SmoothProfile(profile2, profileLen,
            config.slidingWinSize / config.movingStep / 2);
    profile.SmoothProfile(rsProfile1, profileLen,
            config.slidingWinSize / config.movingStep / 2);
    profile.SmoothProfile(rsProfile2, profileLen,
            config.slidingWinSize / config.movingStep / 2);

    LOG_DEBUG2("maxL1Count : "<<maxL1Count);LOG_DEBUG2("maxL2Count : "<<maxL2Count);

    minDist = config.slidingWinSize / config.movingStep + 1;

    getMaxProfilePnt(maxL1Count, profile1);
    getMaxProfilePnt(maxL2Count, profile2);
#ifdef DEBUG
    ranger_debug::dumpArray(profile1, "profile1_n");
    ranger_debug::dumpArray(profile2,"profile2_n");
    ranger_debug::dumpArray(rsProfile1,"rsprofile1_n");
    ranger_debug::dumpArray(rsProfile2,"rsprofile2_n");
#endif
    callPeaks(rsProfile1, rsProfile2, profile1, profile2, config,
            chrom.chromSize, chrom.l1Peaks);
    LOG_DEBUG1("chrom.l1Peaks.size() "<<chrom.l1Peaks.size());
    callPeaks(rsProfile2, rsProfile1, profile2, profile1, config,
            chrom.chromSize, chrom.l2Peaks);
    LOG_DEBUG1("chrom.l2Peaks.size() "<<chrom.l2Peaks.size());

    return chrom.l1Peaks.size();
}
bool CCATPeakFinder::hasLargerNeighbors(size_t minDist, int& tmpStart,
        int& tmpEnd, const vector<size_t>& profile, const peak_t& pk) {
    int tmpStart1, tmpEnd1;
    size_t profileSize = profile.size();
    assert_gt(tmpStart, 0);
    assert_gt(tmpEnd, 0);
    tmpStart1 = tmpStart < minDist ? 0 : tmpStart - minDist;
    tmpEnd1 =
            tmpEnd + minDist + 1 > profileSize ?
                    profileSize - 1 : tmpEnd + minDist;
    LOG_DEBUG5("i, profile[i]: tmpStart1 = "<<tmpStart1);LOG_DEBUG5("i, profile[i]: tmpEnd1 = "<<tmpEnd1);
    //+1 is to match the reverse order in the original code
    assert(tmpEnd+1> -1);
    size_t me = *max_element(profile.begin() + tmpEnd + 1,
            profile.begin() + tmpEnd1 + 1);
    size_t ms = *max_element(profile.begin() + tmpStart1,
            profile.begin() + tmpStart);
    if (me > profile.at(pk.peak) //Must be me >= profile.at()
    || ms >= profile.at(pk.peak)) {
        return true;
    }

    return false;
}
//GetLocalMaxima: find the local maxima in the profile
int CCATPeakFinder::GetLocalMaxima(const vector<size_t>& profile,
        vector<peak_t>& peaks, const int minDist, const size_t minCount) {
    MARK_FUN("CCATPeakFinder::GetLocalMaxima");

    LOG_DEBUG2("minDist "<<minDist);LOG_DEBUG2("minCount "<<minCount);

    int tmpStart, tmpEnd;
    size_t i;
    size_t peakNum;
    bool hasBigdog;
    peakNum = 0;

    assert_geq(peaks.size(), profile.size());
#ifdef DEBUG
    ranger_debug::dumpArray(profile, "profile_n");
#endif
    tmpStart = -1;
    tmpEnd = -1;
    size_t profileSize = profile.size();

    for (i = 1; i < profileSize; i++) {
        LOG_DEBUG5("i, profile[i]: "<<i <<", "<<profile[i]);
        if (profile[i] >= minCount) {
            LOG_DEBUG5("i, profile[i]: profile[i] >= minCount("<<minCount<<")");
            if (profile[i] > profile[i - 1]) {
                LOG_DEBUG5("i, profile[i]: profile[i] > profile[i - 1]("<<profile[i]<<">"<<profile[i - 1]<<")");
                tmpStart = i;
                tmpEnd = i;
                LOG_DEBUG5("i, profile[i]: so that tempStart=tmpEnd=i("<<i<<")");
            }

            if (profile[i] == profile[i - 1]) {
                LOG_DEBUG5("i, profile[i]: profile[i] == profile[i - 1]("<<profile[i]<<"=="<<profile[i - 1]<<")");

                LOG_DEBUG5("i, profile[i]: so that only tmpEnd=i("<<i<<")");
                tmpEnd = i;
            }

            if (profile[i] < profile[i - 1]) {
                LOG_DEBUG5("i, profile[i]: profile[i] < profile[i - 1]("<<profile[i]<<"<"<<profile[i - 1]<<")");
                if ((tmpStart == -1) || (tmpEnd == -1)) {
                    LOG_DEBUG5("i, profile[i]: since (tmpStart == -1) || (tmpEnd == -1), continued.");
                    continue;
                }

                peaks[peakNum].peak = (tmpStart + tmpEnd) / 2;
                LOG_DEBUG5("i, profile[i]: peaks[peakNum].peak =peaks["<<peakNum<<"].peak = "<<(tmpStart + tmpEnd) / 2);

                if (hasLargerNeighbors(minDist, tmpStart, tmpEnd, profile,
                        peaks[peakNum])) {

                    hasBigdog = true;
                    tmpStart = -1;
                    tmpEnd = -1;

                    continue;
                }

                LOG_DEBUG5("Added peak at "<< i);
                peakNum++;
                LOG_DEBUG5("i, profile[i]: peakNum++");
                tmpStart = -1;
                tmpEnd = -1;
            }
        } else {
            LOG_DEBUG5("i, profile[i]: profile[i] < minCount("<<minCount<<")");
            if ((tmpStart == -1) || (tmpEnd == -1)) {
                continue;
            }

            peaks[peakNum].peak = (tmpStart + tmpEnd) / 2;

            if (hasLargerNeighbors(minDist, tmpStart, tmpEnd, profile,
                    peaks[peakNum])) {
                tmpStart = -1;
                tmpEnd = -1;
                continue;
            }

            LOG_DEBUG5("Added peak at "<< i);
            peakNum++;
            tmpStart = -1;
            tmpEnd = -1;
        }
    }

    peaks.resize(peakNum);
    return peakNum;
}
} /* namespace ccat_aux */
