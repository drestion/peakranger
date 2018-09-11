/*
 * fdr_based_thresholder.h
 *
 *  Created on: July 27, 2011
 *      Author: xin
 */

#ifndef FDR_BASED_THRESHOLDER_H_
#define FDR_BASED_THRESHOLDER_H_

#include "region_detector.h"
#include "region_profile/RegionProfile.h"
#include "short_reads/reads.h"
#include "option_parser/cmd_option_parser.h"
#include "result_reporter/result_reporter.h"
#include "utils/logger.h"
#include "utils/assert_helpers.h"
#include "region_profile/profile_summit_detector.h"

#include <vector>
#include <utility>
#include <map>
#include <stdint.h>
#include <iostream>
#include <ostream>

struct thresh {
    // This flag tells whether a stop or start position is currently sought.
    // 0 indicates that a start of a peak is sought and 1 indicates that a
    // stop is sought.
    bool flag;
    // This counts the number of peaks found so far.
    uint32_t count;
    // This stores the most-recently-found stop position.  It is used to
    // determine whether the next peak should be merged with the previous peak
    // or if the new peak should be called a new peak.
    uint32_t stop;

    // At the end of all simulations this field will contain the average
    // count.
    float expectedNumOfPeaks;
};

typedef struct thresh Thresh;

typedef std::map<std::string, std::vector<called_peak> > enriched_regions;
typedef std::vector<uint32_t>::iterator read_itr;
typedef std::pair<read_itr, read_itr> readsregion;
typedef std::vector<std::pair<uint32_t, int8_t> > reads_count_t;
struct window {

    int n_reads;
    // An array holding the count of actual peaks at each threshold.  The
    // index of the array is the threshold - the minimum threshold.
    std::vector<std::vector<called_peak> > peaksAtEachTh;
};

typedef struct window Window;

struct actThresh {
    // This flag tells whether a stop or start position is currently sought.
    // 0 indicates that a start of a peak is sought and 1 indicates that a
    // stop is sought.
    bool peakUnfinished;
    // This stores the most recently found stop position.  It is used to
    // determine whether the next peak should be merged with the previous peak
    // or if it is a new peak.
    uint32_t endOfLastPeak;
    //    // This is a pointer to the window indicated by peaks in this threshold.
//    // It is used to count a peak in each window it overlaps.
//    uint32_t currentWindow;

    // An array of the peaks at this threshold.
    std::vector<called_peak> peaks;
    actThresh() {
        peaks.clear();
        peakUnfinished = false;
        endOfLastPeak = 0;
    }

};

typedef struct actThresh ActThresh;

struct _SGR {
    std::vector<uint32_t> poss;
    std::vector<int32_t> heights;
};

typedef struct _SGR SGR;

/*
 * Merge reads on both strands into the pos strands
 *
 * detect regions above the cut-off as peaks.
 *
 * Good for high-quality sharp peaks, not good
 *
 * for histone marks.
 */
class fdr_based_thresholder: public region_detector {
public:
    fdr_based_thresholder() :
            region_detector() {

    }
    ~fdr_based_thresholder(){}
    void detectSummits(Reads& treatment_reads, Reads& control_reads,
            cmd_option_parser& option);

    void detectSummits(Reads& treatment_reads, Reads& control_reads,
            cmd_option_parser& option, std::ostream& os);

    void detectSummits(Reads& treads, Reads& creads, uint32_t no_of_thread,
            double p_val_cutoff, double delta, uint32_t threshold,
            uint32_t mergedistance, uint32_t binlength, uint32_t bandwidth,
            uint32_t readextension, std::ostream& os, bool verbose);
     void export_results(result_reporter& reporter, std::ostream& om){std::cerr<<"Not implemented yet\n";}
protected:

     uint32_t _calulate_number_of_bins(uint32_t p_chr_length, uint32_t n_chr_length, uint32_t binlength) const;
     uint32_t _calculate_pos_read_loci_ret_start(uint32_t start, uint32_t readlength, uint32_t extension ) const;
     uint32_t _calculate_pos_read_loci_ret_end(uint32_t start, uint32_t readlength, uint32_t extension ) const;
     uint32_t _calculate_neg_read_loci_ret_start(uint32_t start, uint32_t readlength, uint32_t extension ) const;
     uint32_t _calculate_last_pos(int32_t pos_read, uint32_t neg_read,
    			uint32_t readlength, uint32_t extension, uint32_t pchrlength, uint32_t nchrlength) const;
         uint32_t _calculate_neg_read_loci_ret_end(uint32_t start, uint32_t readlength, uint32_t extension ) const;


    void _setup(uint32_t no_of_threads, uint32_t binlength, uint16_t threshold,
            uint32_t mergedistance, uint32_t readlength, uint32_t readextension,
            uint32_t bandwidth, double pval, double delta) {

        _nThreads = no_of_threads;
        _binlength = binlength;
        _threshold = threshold;
        _mergedistance = mergedistance;
        _readlength = readlength;
        _readextension = readextension;
        _bandwidth = bandwidth;
        _pvalThreshold = pval;
        _delta = delta;

    }

    void _normalize_reads(Reads& treads, Reads& creads);
    inline double _normalize_reads_onchr(Reads& treads, Reads& creads,
            std::string& chr, uint32_t binlength);

    void _processChr(std::ostream& os, bool print_stream);

    void _get_random_reads(uint32_t start, uint32_t end, uint32_t no_of_reads,
            std::vector<uint32_t>& result) const;

    void getRndPeaksAtEachThreshold(uint32_t start, uint32_t end,
            uint32_t no_of_reads, uint32_t extension,
            std::vector<Thresh>& rndPeaks) const;

    void walk(reads_count_t& reads_count, std::vector<Thresh>& rndPeaks) const;

    uint32_t findThresh(std::vector<Thresh>& rndPeaks,
            std::vector<ActThresh>& realPeaks, double& fdr) const;
    void process(std::vector<ActThresh>& rndPeaks, const uint32_t height,
            const uint32_t pos) const;

    void findPeaksAtEachThreshold(SGR& sgr, std::vector<ActThresh>& realPeaks,
            uint32_t binend) const;

    virtual void get_profile_of_reads(uint32_t extension, uint32_t readlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, SGR & result) const;

protected:
    struct workerPara {
        workerPara():_verbose(false) {
        }

        std::vector<std::string> _chrsleft;
        std::map<std::string, readsregion> _all_pos_treads;
        std::map<std::string, readsregion> _all_neg_treads;
        std::map<std::string, readsregion> _all_pos_creads;
        std::map<std::string, readsregion> _all_neg_creads;
        bool _verbose;
        std::map<std::string, readsregion> getAllNegCreads() const;
        std::map<std::string, readsregion> getAllNegTreads() const;
        std::map<std::string, readsregion> getAllPosCreads() const;
        std::map<std::string, readsregion> getAllPosTreads() const;
        std::vector<std::string> getChrsleft() const;
        bool isVerbose() const;
        void setAllNegCreads(
                std::map<std::string, readsregion> _all_neg_creads);
        void setAllNegTreads(
                std::map<std::string, readsregion> _all_neg_treads);
        void setAllPosCreads(
                std::map<std::string, readsregion> _all_pos_creads);
        void setAllPosTreads(
                std::map<std::string, readsregion> _all_pos_treads);
        void setChrsleft(std::vector<std::string> _chrsleft);
        void setVerbose(bool _verbose);
    };

    workerPara _workerPara;
protected:
    uint32_t _nThreads;
    uint32_t _binlength;
    uint16_t _threshold;
    uint32_t _mergedistance;
    uint32_t _readlength;
    uint32_t _readextension;
    uint32_t _bandwidth;
    double _pvalThreshold;
    double _delta;
    std::map<std::string, double> _norm_factors;

}
;

#endif /* FDR_BASED_THRESHOLDER_H_ */
