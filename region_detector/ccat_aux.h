/*
 * ccat_aux.h
 *
 *  Created on: Mar 27, 2012
 *      Author: xin
 */

#ifndef CCAT_AUX_H_
#define CCAT_AUX_H_

#include "ccat_peak.h"
#include "chrt.h"
#include "ccat_config.h"
#include "ccat_bint.h"
#include "short_reads/reads.h"
#include <algorithm>
#define MIN_COUNT 10000
#define BIN_SIZE 1000
#define MAX_ITERATION 20
#define TRANSITION_SIZE 30
#define Q_VALUE_STEP_NUM 1000

namespace ccat_aux {
void initializeSrandUsingCurrentTime();
void LoadData(std::vector<ccat_aux::chr_t> & chroms, Reads & treads,
        Reads & creads, size_t & chromNum);
void SortAndDedup(std::vector<ccat_aux::chr_t> & chroms);
void load_histone_config(ccat_aux::ccat_config_t & config);
void config_validation(ccat_aux::ccat_config_t & config);

double ComputeFoldChange(size_t count1, size_t count2, size_t tagCount,
        size_t binCount, double smoothingFactor);

//GenRegionProfile: generate profile of significant regions.
void GenRegionProfile(std::vector<size_t> & l1Profile,
        std::vector<size_t> & l2Profile, std::vector<size_t> & region,
        double l1Ratio, double l2Ratio, size_t slidingWinSize,
        size_t movingStep, size_t bootstrapPass, double minScore,
        size_t l1MaxCount, size_t l2MaxCount, size_t tagCount, size_t binCount,
        double smoothingFactor, std::vector<double>& lookUpTable,
        std::vector<int>& flag);

double BootstrapFoldChange(size_t l1Count, size_t l2Count, double l1Ratio,
        double l2Ratio, size_t bootstrapPass, size_t tagCount, size_t binCount,
        double smoothingFactor);

double ComputeSmoothingParameter(std::vector<peak_t>& peaks, size_t& binCount,
        size_t& tagCount);

//fit a negative binomial distribution
int FitNegBinomDist(std::vector<int>& hist, double& m, double& k, double k_min,
        double k_max);

//Compute log-gamma function
double gammln(double xx);

template<typename Container>
size_t SortAndGetUniqueTags(Container& tags) {

    std::sort(tags.begin(), tags.end());

    tags.resize(std::unique(tags.begin(), tags.end()) - tags.begin());

    return tags.size();
}

class mm {
public:

    static std::ostream& map(std::ostream& os, const peak_t& pk) {
        os << pk.foldChange;
        return os;
    }
    static std::string map(const peak_t& pk) {
        std::stringstream ss;
        ss << pk.foldChange;
        return ss.str();
    }
};

class bin_t_l1 {
public:
    static int& map(bin_t& pk) {

        return pk.l1Counts;
    }
};

class bin_t_l2 {
public:
    static int& map(bin_t& pk) {
        return pk.l2Counts;
    }
};

}

#endif /* CCAT_AUX_H_ */
