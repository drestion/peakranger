/*
 * CCATfdrAnalyzer.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef CCATFDRANALYZER_H_
#define CCATFDRANALYZER_H_

#include "common/stl_header.h"
#include "ccat_peak.h"
#include "ccat_profile.h"
#include "chrt.h"
#include "ccat_config.h"
namespace ccat_aux {

class CCATfdrAnalyzer {
public:
    CCATfdrAnalyzer();
    virtual ~CCATfdrAnalyzer();

    //SignificanceAnalysis: perform significance analysis, will be called by the main routine of CCAT
    int SignificanceAnalysis(std::vector<ccat_aux::chr_t>&chroms, size_t chromNum, double l1Ratio, double l2Ratio,
            size_t maxL1Count, size_t maxL2Count, const ccat_aux::ccat_config_t& config);

private:

    std::vector<double> q;
    std::vector<double> value;
    std::vector<double> lookUpTable;
    std::vector<int> flag;
    size_t row;
    size_t column;
    size_t tagCount;
    size_t binCount;
    double smoothingFactor;
    ccat_profile_t profile;
    const static size_t QVALUESTEP;

    //ComputeThreshold: compute the local FDR, return the FDR at the cut-off threshold;
    double ComputeLocalFDR(std::vector<ccat_aux::peak_t>& l1Peaks, std::vector<ccat_aux::peak_t>& l2Peaks,
            std::vector<double> &q, std::vector<double>&value, const ccat_aux::ccat_config_t& config, size_t tagCount,
            size_t binCount);
    //PostProcessing: post-processing of the peaks: boostrapping for fold-change calculation, region identification, peak refinement
    int PostProcessing(std::vector<ccat_aux::chr_t>& chroms, int chromNum, double l1Ratio, double l2Ratio,
            int maxL1Count, int maxL2Count, const ccat_aux::ccat_config_t& config);

    //ProcessOneChrom: process one chromosome
    int ProcessOneChrom(ccat_aux::chr_t& chrom, double l1Ratio, double l2Ratio, const ccat_aux::ccat_config_t& config);

};

} /* namespace ccat_aux */
#endif /* CCATFDRANALYZER_H_ */
