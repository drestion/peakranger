/*
 * CCATNoiseRate.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef CCATNOISERATE_H_
#define CCATNOISERATE_H_

#include "chrt.h"
#include "ccat_bint.h"

namespace ccat_aux {

class CCATNoiseRate {
public:
    CCATNoiseRate();
    virtual ~CCATNoiseRate();

    //ComputeNoiseRate: compute the noise rate
    double ComputeNoiseRate(const std::vector<chr_t>& chroms, size_t chromNum,
            double& l1Ratio, double& l2Ratio) const;
    void setFragmentSize(const size_t& fs);

private:

    size_t fragSize;

    void ComputeSampleRatio(double noiseRate, int l1TagNum, int l2TagNum,
            double & l1Ratio, double & l2Ratio) const;
    //ComputeNoiseRateInOneIteration: compute the noise rate by one iteration
    void ComputeNoiseRateInOneIteration(const std::vector<chr_t> & chroms,
            size_t chromNum, double & noiseRate) const;
    //CountFragmentsInOneChrom: sample portion of the tags by ratio, and add fragment counts to the bins
    void CountFragmentsInOneChrom(const chr_t & chrom,
            std::vector<bin_t> & bins, int binNum, double l1Ratio,
            double l2Ratio, size_t fragmentSize) const;

    void AllocBinMem(const std::vector<chr_t> & chroms, int chromNum,
            std::vector<std::vector<bin_t> >& bins,
            std::vector<size_t> & binNum) const;

    const static size_t BINSIZE;
    const static size_t MAXITERATION;
    const static size_t MINCOUNT;

};

} /* namespace ccat_aux */
#endif /* CCATNOISERATE_H_ */
