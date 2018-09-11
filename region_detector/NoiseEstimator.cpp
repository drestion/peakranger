/*
 * NoiseEstimator.cpp
 *
 *  Created on: Jun 12, 2012
 *      Author: xfeng
 */

#include "region_detector/NoiseEstimator.h"
#include "region_detector/ccat_aux.h"

using namespace ccat_aux;
using namespace std;
namespace reads {

NoiseEstimator::NoiseEstimator() :
        mNR(), mFragSize(200), mSeed(123456) {

}

NoiseEstimator::~NoiseEstimator() {

}

double NoiseEstimator::estimate(Reads& treads, Reads& creads) {

    vector<chr_t> chroms;
    size_t chromNum;
    double l1Ratio = 1;
    double l2Ratio = 1;

    LoadData(chroms, treads, creads, chromNum);
    SortAndDedup(chroms);
    initializeSrandUsingCurrentTime();
    mNR.setFragmentSize(mFragSize);

    return mNR.ComputeNoiseRate(chroms, chromNum, l1Ratio, l2Ratio);
}

} /* namespace reads */
