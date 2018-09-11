/*
 * NoiseEstimator.h
 *
 *  Created on: Jun 12, 2012
 *      Author: xfeng
 */

#ifndef NOISEESTIMATOR_H_
#define NOISEESTIMATOR_H_
#include "region_detector/CCATNoiseRate.h"
#include "short_reads/reads.h"
namespace reads {

class NoiseEstimator {
public:
    NoiseEstimator();
    virtual ~NoiseEstimator();

    double estimate(Reads& treads, Reads& creads);

    uint32_t getFragSize() const {
        return mFragSize;
    }

    void setFragSize(uint32_t fragSize) {
        mFragSize = fragSize;
    }

    int32_t getSeed() const {
        return mSeed;
    }

    void setSeed(int32_t seed) {
        mSeed = seed;
    }

private:
    ccat_aux::CCATNoiseRate mNR;
    uint32_t mFragSize;
    int32_t mSeed;
};

} /* namespace reads */
#endif /* NOISEESTIMATOR_H_ */
