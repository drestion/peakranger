/*
 * ccat.h
 * A wrapper for CCAT.
 *
 * Han Xu et al.
 * A signal noise model for significance analysis of ChIP-seq with negative control
 * Bioinformatics (2010) 26 (9): 1199-1204.
 *
 *  Created on: Dec 14, 2011
 *      Author: xfeng
 */

#ifndef CCAT_H_
#define CCAT_H_
#include "region_detector.h"
#include "ccat_aux.h"
#include "ccat_config.h"
#include "ccat_profile.h"
#include "CCATNoiseRate.h"
#include "CCATfdrAnalyzer.h"
#include "CCATPeakFinder.h"
#include <iostream>
class ccat: public region_detector {
public:
    ccat();
    virtual ~ccat();
    virtual void detectSummits(Reads& treatment_reads, Reads& control_reads,
            cmd_option_parser& option);

    virtual void detectSummits(Reads& treatment_reads, Reads& control_reads,
            cmd_option_parser& option, std::ostream& os) {
    }

private:
    int UploadPeaks(std::vector<ccat_aux::chr_t>&chroms, size_t chromNum,
            const char *projectName, const ccat_aux::ccat_config_t& config);

    void insertPeak(const std::string& chr, called_peak& pk);

    void cmain(Reads& treads, Reads& creads, cmd_option_parser& option);

    void LoadData(std::vector<ccat_aux::chr_t> & chroms, Reads & treads,
            Reads & creads, size_t & chromNum);
    void SortAndDedup(std::vector<ccat_aux::chr_t> & chroms);

    std::vector<double> lookUpTable;
    std::vector<int> flag;
    std::vector<double> q;
    std::vector<double> value;

    ccat_aux::ccat_config_t conf;
    ccat_aux::ccat_profile_t prof;
    ccat_aux::CCATNoiseRate noise;
    ccat_aux::CCATfdrAnalyzer fdr;
    ccat_aux::CCATPeakFinder pkFinder;
};

#endif /* CCAT_H_ */
