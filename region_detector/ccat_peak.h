/*
 * ccat_peak.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef CCAT_PEAK_H_
#define CCAT_PEAK_H_

#include "common/stl_header.h"

namespace ccat_aux {

class peak_t{
    friend std::ostream& operator<<(std::ostream& os,
                                    peak_t& pk);
    public:
    peak_t();
    virtual ~peak_t();
    size_t chromIndex; //Index of chromosome
    int start; //start of the peak region
    int end; //end
    int peak; //peak location
    int l1Count; //read counts in L1
    int l2Count; //read counts in L2
    int reSampledL1Count; //resampled read counts in L1
    int reSampledL2Count; //resampled read counts in L2
    double foldChange; //fold change
    double qValue; //q value
    int isSignificant; //1 if the peak is signficant, otherwise 0


};
bool lessFoldChange(const peak_t& lhs,
                    const peak_t& rhs);
bool greaterFoldChange(const peak_t& lhs,
                       const peak_t& rhs);
} /* namespace ccat_aux */

#endif /* CCAT_PEAK_H_ */

