/*
 * fdr_based_thresholder
 *
 *  Created on: July 27, 2011
 *      Author: xin
 */

#include "fdr_based_thresholder.h"
#include "region_profile/RegionProfile.h"
#include "math/distributions.h"
#include "region_profile/arrayThresholder.h"
#include "region_profile/points_profile.h"
#include "utils/exceptions.h"
#include "region_profile/profile_smoother.h"
#include "region_profile/peakseq_profile_thresholder.h"
#include "region_profile/peakseq_profile.h"
#include "option_parser/cmd_option_parser.h"
#include "utils/timer.h"
#include "utils/assert_helpers.h"
#include "short_reads/point_reads.h"
#include "wiggle/wigbuilder.h"
#include "wiggle/wig.h"
#include "short_reads/reads_aux.h"

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include <math.h>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>
#include <sstream>

#define MERGE_DISTANCE 200
#define SIM_ROUND      3
#define MAX_THRESHOLD  100
#define MIN_THRESHOLD  2
#define NUM_THRESHOLD  MAX_THRESHOLD - MIN_THRESHOLD + 1 //99
#define MIN_FDR  0.001
#define MIN_PEAKSIZE   2
#define MAX_PEAKSIZE   100000
#define foreach BOOST_FOREACH
using namespace std;
typedef vector<pair<uint32_t, int8_t> > reads_count_t;

namespace {
  boost::mutex _chr_mt;
  boost::mutex _result_mt;

  template<typename T>
    inline void summits_2_vector(vector<T>& vec, T offset, vector<T>& result) {
      size_t size = 0;
      while (size < vec.size()) {
        result.push_back(vec[size++] + offset);
      }
    }

  /*
   * parameter struct shared by all workers
   */

  uint32_t trans(uint32_t read, uint32_t readlength, uint32_t extension) {
    uint32_t sum = read + readlength;
    if (sum > extension) {
      return sum - extension;
    } else
      return 0;
  }


  bool sort_comparator_int8(pair<uint32_t, int8_t> p1,
      pair<uint32_t, int8_t> p2) {
    return p1.first < p2.first;
  }

  void get_reads_count(vector<uint32_t> reads, uint32_t extension,
      reads_count_t& reads_count) {
    LOG_DEBUG1("fdr_based_thresholder:get_reads_count");
    uint32_t a = 0;
    uint32_t b = 0;
    uint32_t read = 0;

    typedef vector<pair<uint32_t, int8_t> > reads_count_t;
    vector<uint32_t>::iterator readsStart, readsEnd;
    readsStart = reads.begin();
    readsEnd = reads.end();

    while (readsStart != readsEnd) {
      read = *readsStart++;
      a = read;
      b = read + extension;
      reads_count.push_back(pair<uint32_t, int8_t>(a, 1));
      reads_count.push_back(pair<uint32_t, int8_t>(b, -1));
    }

    //sort based on their locations;
    sort(reads_count.begin(), reads_count.end(), sort_comparator_int8);

    LOG_DEBUG1("QUIT:fdr_based_thresholder:get_reads_count");
  }

}/* Namespace */

/*
 * Generates expected number of peaks at each threshold
 * and store them in allThresholds
 */
void fdr_based_thresholder::getRndPeaksAtEachThreshold(uint32_t start,
    uint32_t end, uint32_t no_of_reads, uint32_t extension,
    vector<Thresh>& rndPeaks) const {
  LOG_DEBUG2("fdr_based_thresholder::getRndPeaksAtEachThreshold");
  assert_gt(end, start);
  for (size_t i = 0; i < SIM_ROUND; i++) {
    // Get the sorted array of simulated positions.
    vector<uint32_t> reads;

    _get_random_reads(start, end, no_of_reads, reads);
    LOG_DEBUG3("got "<<reads.size()<<" rand reads");
    reads_count_t rand_reads_count;

    get_reads_count(reads, extension, rand_reads_count);
    LOG_DEBUG3("got "<<rand_reads_count.size()<<" rand reads counts");
    // Walk down the positions beginning at the lowest position and ending at
    // the highest position keeping track of height and cataloguing when each
    // threshold is exceeded.
    walk(rand_reads_count, rndPeaks);

    // For each thresh structure, add the count to the numerator of the
    // average field.
    assert_lt(NUM_THRESHOLD-1, rndPeaks.size())
      for (size_t i = 0; i < NUM_THRESHOLD; i++) {
        rndPeaks[i].expectedNumOfPeaks += (1.0 * rndPeaks[i].count)
          / SIM_ROUND;
      }
  }
#ifdef USE_LOGGING
  assert_lt(70,rndPeaks.size())
    for(size_t i = 0; i < 70; i++) {
      LOG_DEBUG2("Expect " << rndPeaks[i].expectedNumOfPeaks<<" peaks at threshold " <<i);
    }
#endif
  LOG_DEBUG2("QUIT:fdr_based_thresholder::getRndPeaksAtEachThreshold");
}

/*
 * Adapted from original PeakSeq C source codes.
 */
void fdr_based_thresholder::walk(reads_count_t& reads_count,
    vector<Thresh>& rndPeaks) const {

  if (rndPeaks.size() < NUM_THRESHOLD) {
    rndPeaks.resize(NUM_THRESHOLD);
  }
  reads_count_t::iterator cntStart, cntMed, cntEnd;
  cntStart = reads_count.begin();
  cntEnd = reads_count.end();
  cntMed = cntStart;
  int32_t height = 0; // Initialize the height.

  assert_lt(MAX_THRESHOLD-MIN_THRESHOLD, rndPeaks.size())
    // This loop runs once for each element in the array.
    for (size_t i = 0; i < reads_count.size(); i++) {
      // Increment the height by the score at this position.  Thus the
      // height field will always contain the number of reads at this
      // position.
      height += (int32_t) reads_count[i].second;

      // This loop continues while the position is constant.  In this way,
      // position structures with the same location are treated at the same
      // time.

      while (i + 1 < reads_count.size()) {
        if (reads_count[i].first == reads_count[i + 1].first) {
          i++;
          height += (int32_t) reads_count[i].second;
        } else {
          break;
        }
      }
      // This loop runs once for each threshold.
      for (uint32_t thresh = MIN_THRESHOLD; thresh <= MAX_THRESHOLD;
          thresh++) {
        // Get a pointer to the current threshold's structure.
        Thresh& t = rndPeaks[thresh - MIN_THRESHOLD];

        // If the height is above the threshold and we are looking for a
        // starting position of a peak, count the peak.
        if (height > 0) {
          if ((height >= (int32_t) thresh) && (t.flag == false)) {
            // If this peak is further than the maximum gap from the most
            // recently found peak, then this peak is a new peak and must
            // be counted.
            if (reads_count[i].first - t.stop > MERGE_DISTANCE) {
              LOG_DEBUG4("In fdr_based_thresholder::walk, got 1 peaks at threshold "
                  <<thresh <<" for the " <<i<<"th reads_count");
              t.count++;
            }

            // Whether this is a new peak or not, set the flag to 1 to
            // look for the end of this peak.
            t.flag = true;
          }
          // If the height is below the threshold and the end of a peak is
          // sought, then this is the end of the peak.
          else if ((t.flag == true) && (height < (int32_t) thresh)) {
            t.stop = cntStart->first;
            t.flag = false;
          }
        }
      }
    }

#ifdef USE_LOGGING
  for(size_t i = 0; i < 5; i++) {
    LOG_DEBUG2("Finally, In fdr_based_thresholder::walk, got " << rndPeaks[i].count <<" peaks at threshold " <<i);
  }
#endif
}

void fdr_based_thresholder::_get_random_reads(uint32_t start, uint32_t end,
    uint32_t no_of_reads, vector<uint32_t>& reads) const {
  assert_gt(end, start)
    // distributions::random_int_std(start, end, no_of_reads, reads);
    distributions::random_int_boost(start,end,no_of_reads,reads);

  sort(reads.begin(), reads.end());
}
/*
 * Adapted from original PeakSeq C source codes.
 */
void fdr_based_thresholder::findPeaksAtEachThreshold(SGR& sgr,
    vector<ActThresh>& realPeaks, uint32_t binend) const {
  //The SGR is usually a region's SGR
  LOG_DEBUG1("fdr_based_thresholder::findPeaksAtEachThreshold");
  uint32_t pos;
  int32_t height;

  if (sgr.poss.size() < 1)
    return;

  if (realPeaks.size() < NUM_THRESHOLD) {
    realPeaks.resize(NUM_THRESHOLD);
    realPeaks.assign(NUM_THRESHOLD, ActThresh());
  }


  //Iterate all signal profiles and find peaks
  //using all defiend thresholds, e.g 0-1000 
  for (size_t i = 0; i < sgr.poss.size(); i++) {
    pos = sgr.poss[i];
    height = sgr.heights[i];
    process(realPeaks, height, pos);
  }

  // Store the stop position of the last peak as long as there was at least
  // one peak.

  //This loop works to finish unfinished peaks.

  for (int i = 0; i < NUM_THRESHOLD; i++) {
    // If there was at least 1 peak at this threshold, the following
    // comparisons will be safe.
    if (realPeaks[i].peaks.size() > 0) {
      // If the stop of this peak was recorded in the at.second field,
      // use at.second.
      // This happens when no reads are below the threshold for the peak.
      // TODO: is this normal?
      if (realPeaks[i].endOfLastPeak > realPeaks[i].peaks.back().first) {
        realPeaks[i].peaks.back().second = realPeaks[i].endOfLastPeak;
      }

      // Otherwise the peak continued like a flat plateau,
      // use the boundary as the end of the peak
      // This should be normal, although rare.
      // TODO: merge these peaks later? two adjacent bins will have peaks close to each
      else {
        //TODO:A better boundary should be the last read's start or last sgr
        //TODO:May be should just control this from the caller function
        realPeaks[i].peaks.back().second = binend;
      }
    }
  }
#ifdef USE_LOGGING
  for(size_t i = 0; i < 5; i++) {
    LOG_DEBUG2("Found " << realPeaks[i].peaks.size()<<" real peaks at threshold " <<i);
  }
#endif
  LOG_DEBUG1("QUITED:fdr_based_thresholder::findPeaksAtEachThreshold");
}
/*
 * Adapted from original PeakSeq C source codes.
 */
/*
 * realPeaks.size() > 1
 */
void fdr_based_thresholder::process(vector<ActThresh>& realPeaks,
    const uint32_t height, const uint32_t pos) const {

  if (realPeaks.size() < MAX_THRESHOLD - MIN_THRESHOLD + 1) {
    realPeaks.resize(MAX_THRESHOLD - MIN_THRESHOLD + 1);
  }

  LOG_INFO("fdr_based_thresholder::process");
  LOG_DEBUG3("got height: "<<height <<" at position "<<pos);

  // This loop runs once for each threshold.
  for (uint32_t th = MIN_THRESHOLD; th <= MAX_THRESHOLD; th++) {
    // Get the threshold index.
    uint32_t actInd = th - MIN_THRESHOLD;

    // Get the current threshold's structure.
    ActThresh& at = realPeaks[actInd];

    // This block is entered if the height is above the threshold and we
    // are looking for the starting position of a peak.
    if ((height >= th) && !(at.peakUnfinished)) {
      // If this peak is further than the maximum gap from the most
      // recently found peak, then this peak is a new peak and must be
      // counted.
      if (pos - at.endOfLastPeak > MERGE_DISTANCE) {

        if (at.peaks.size() > 0) {
          //Mark the boundary of the previously sought peak
          at.peaks.back().second = at.endOfLastPeak;
        }

        LOG_DEBUG4("got the begin of a new peak with height: "
            <<height <<" at position "<<pos);

        // Insert the new peak into the array.

        called_peak tmp;
        tmp.first = pos;
        tmp.second = pos;
        tmp.p = -1;
        at.peaks.push_back(tmp);
      }
      else {
        //The peak is merged
      }
      // Whether this is a new peak or not, set the flag to 1 to look
      // for the end of this peak.
      at.peakUnfinished = true;
    }
    // If the height is below the threshold and the end of a peak is
    // being sought, this is the end of the peak.
    // It is possible that the height reached the boundary of the bin
    // leaving an incomplete peak.
    else if ((height < th) && (at.peakUnfinished)) {
      // Set the stop field to this position and set the flag to 0.  If
      // the next peak begins more than max_gap away from this peak,
      // then this peak's stop position will be set to the value of the
      // at.stop field.
      //TODO: -1 or not?
      //            at.endOfLastPeak = pos - 1;
      if (pos > 1) {
        at.endOfLastPeak = pos - 1;
      } else {
        at.endOfLastPeak = 0;
      }
      at.peakUnfinished = false;
    }
    //else if ((height >= th) && (at.peakUnfinished)) {
    //These unfinished peaks will be finished in the mother function
    //}
    //else if ((height < th) && !(at.peakUnfinished)) {
    //Otherwise, the height is just skipped
    //}
  }
  LOG_DEBUG2("QUITED:fdr_based_thresholder::process");
}


uint32_t fdr_based_thresholder::_calculate_last_pos(int32_t pos_read,
		uint32_t neg_read, uint32_t readlength, uint32_t extension,
		uint32_t pchrlength, uint32_t nchrlength) const {
	 if(pchrlength >= nchrlength){
		             	  return _calculate_pos_read_loci_ret_end(pos_read,readlength,extension);
		               }
		               else{
		            	   return _calculate_neg_read_loci_ret_end(neg_read,readlength,extension);
		               }
}

// Input reads do not have to be sorted
void fdr_based_thresholder::get_profile_of_reads(uint32_t extension,
    uint32_t readlength, vector<uint32_t>::iterator readsStart,
    vector<uint32_t>::iterator readsEnd,
    vector<uint32_t>::iterator nreadsStart,
    vector<uint32_t>::iterator nreadsEnd, SGR & result) const {
  uint32_t a = 0;
  uint32_t b = 0;
  uint32_t read = 0;

  typedef vector<pair<uint32_t, int8_t> > reads_count_t;
  reads_count_t reads_count;

  LOG_DEBUG1("fdr_based_thresholder::get_profile_of_reads");
  /*
   * Positive reads
   */
  while (readsStart != readsEnd) {
    read = *readsStart++;
    a = _calculate_pos_read_loci_ret_start(read,readlength,extension);
    b = _calculate_pos_read_loci_ret_end(read,readlength,extension);
    reads_count.push_back(pair<uint32_t, int8_t>(a, 1));
    reads_count.push_back(pair<uint32_t, int8_t>(b, -1));
  }
  /*
   * Negative reads
   */
  while (nreadsStart != nreadsEnd) {
    read = *nreadsStart++;
    if (read + readlength < extension)
      continue;

    a = _calculate_neg_read_loci_ret_start(read,readlength,extension);
    b = _calculate_neg_read_loci_ret_end(read,readlength,extension);
    reads_count.push_back(pair<uint32_t, int8_t>(a, 1));
    reads_count.push_back(pair<uint32_t, int8_t>(b, -1));
  }
  if (reads_count.size() == 0) {
    return;
  }
  /*
   * Sort based on their locations;
   */

  sort(reads_count.begin(), reads_count.end(), sort_comparator_int8);
#ifdef USE_LOGGING
  for(int i = 0; i < reads_count.size(); i++)
    cout <<reads_count[i].first<<"\t"<<static_cast<int>(reads_count[i].second)<<endl;
#endif
  int32_t height = (int32_t) reads_count.front().second;
  uint32_t last_pos = reads_count.front().first;
  // This loop runs over all pos structures in this array.
  for (size_t i = 1; i < reads_count.size(); i++) {

    if (reads_count[i].first == last_pos) {
      height += (int32_t) reads_count[i].second;
      continue;
    }
    if (last_pos > 0 && height > 0) {
      result.heights.push_back(height);
      result.poss.push_back(last_pos);
    }
    height += reads_count[i].second;
    last_pos = reads_count[i].first;

  }
  if (last_pos > 0 && height > 0) {
    result.heights.push_back(height);
    result.poss.push_back(last_pos);
  }

  LOG_DEBUG1("FINISHED building sgr, total profile points: "<<result.poss.size());
  LOG_DEBUG1("QUITED:fdr_based_thresholder::get_profile_of_reads");
}
/*
 * Adapted from original PeakSeq C source codes.
 */
/*
 * rndPeaks.size() > 0;
 * realpeaks.size() > 0;
 *
 *
 */
uint32_t fdr_based_thresholder::findThresh(vector<Thresh>& rndPeaks,
    vector<ActThresh>& realPeaks, double& fdr) const {
  LOG_DEBUG1("fdr_based_thresholder::findThresh");
  // Check if the minimum threshold qualifies.
  /// Calculate the FDR at the minimum threshold.
#ifdef USE_LOGGING
  for(size_t i = 0; i < 15; i++) {
    LOG_DEBUG2("realPeaks: " << realPeaks[i].peaks.size()<<"\trndPeaks: " <<rndPeaks[i].expectedNumOfPeaks);
  }
#endif
  if (realPeaks[0].peaks.size() < 1){
    //If no peaks were called even at zero...
    //Avoid overflow
    //This is usually very rare
    return MIN_THRESHOLD ;
  }
  fdr = rndPeaks[0].expectedNumOfPeaks / realPeaks[0].peaks.size();
  // If the FDR is below the minimum threshold and the average count is
  // decreasing, the minimum threshold qualifies.
  // This is usually extremely rare giving the capacity of modern sequencers
  // and fdr is now usually larger than 1
  if ((fdr <= MIN_FDR)
      && (rndPeaks[1].expectedNumOfPeaks <= rndPeaks[0].expectedNumOfPeaks)) {
    return MIN_THRESHOLD;
  }

  // Initialize the threshold index to the index of the second-lowest
  // threshold.
  size_t th = 1;

  // Advance until the counts begin decreasing.
  while (th < NUM_THRESHOLD) {
    if (rndPeaks[th].expectedNumOfPeaks
        >= rndPeaks[th - 1].expectedNumOfPeaks) {
      // Usually the line is straight down so this iterates only once
      // And thus it is really hard to reach here
      if (rndPeaks[th].expectedNumOfPeaks != 0) {
        th++;
      } else {
        //TODO: how to reach this point?
        return th + MIN_THRESHOLD;
      }
    } else {
      break;
    }
  }LOG_DEBUG1("fdr_based_thresholder::found descending threshold at "<< th);
  // Once the counts are decreasing, find the first value that qualifies as
  // below the minimum FDR.
  for (; th < NUM_THRESHOLD; th++) {
    // Correct for the case that the actual count (the denominator) is 0.
    if (realPeaks[th].peaks.size() == 0) {
      return th + MIN_THRESHOLD;
    }
    // Otherwise, calculate the FDR at this threshold.
    fdr = rndPeaks[th].expectedNumOfPeaks / realPeaks[th].peaks.size();

    // If the FDR is below the required FDR, then this is the threshold.
    if (fdr <= MIN_FDR){
      return (th + MIN_THRESHOLD); }
  }

  // No threshold was found < MAX_THRESHOLD
  return MAX_THRESHOLD;
}

void fdr_based_thresholder::_processChr(ostream& os, bool print_stream) {
  typedef pair<vector<uint16_t>::iterator, vector<uint16_t>::iterator> profile_region;
  vector<profile_pos_t> profile;
  vector<uint16_t> null_profile(0, 0);

  uint32_t a, b, no_of_reads;
  double p, norm_factor, fdr;
  string chr;

  map<string, readsregion> __all_pos_treads;
  map<string, readsregion> __all_neg_treads;
  map<string, readsregion> __all_pos_creads;
  map<string, readsregion> __all_neg_creads;
  bool __verbose;
  bool _profilenotgood = false;
  {
    boost::mutex::scoped_lock l(_chr_mt);
    if (_workerPara._chrsleft.empty()) {
      //no more chrs to process
      return;
    }
    chr = _workerPara._chrsleft.back();
    //        cout <<"norm_factor:"<<chr<<":"<<norm_factor<<endl;
    _workerPara._chrsleft.pop_back();
    norm_factor = _norm_factors[chr];

    __all_pos_treads = _workerPara._all_pos_treads;
    __all_neg_treads = _workerPara._all_neg_treads;
    __all_pos_creads = _workerPara._all_pos_creads;
    __all_neg_creads = _workerPara._all_neg_creads;
    __verbose = _workerPara._verbose;

  }
  uint32_t __binlength = _binlength;
  uint32_t __readlength = _readlength;
  uint32_t __readextension = _readextension;
  uint32_t __bandwidth = _bandwidth;
  double __pvalThreshold = _pvalThreshold;
  double __delta = _delta;

  //index for reads in the chromosome(p), bin(pp) and the peak (ppp)
  vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
    preadsend, ppreadsend, pppreadsend, pcreadsstart, ppcreadsstart,
    pppcreadsstart, pcreadsend, ppcreadsend, pppcreadsend;

  vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
    npreadsend, nppreadsend, npppreadsend, npcreadsstart,
    nppcreadsstart, npppcreadsstart, npcreadsend, nppcreadsend,
    npppcreadsend;

  //Loop until all chrs have been processed
  while (true) {

    //get reads in this chr
    try {
      preadsstart = __all_pos_treads[chr].first;
      preadsend = __all_pos_treads[chr].second;
      pcreadsstart = __all_pos_creads[chr].first;
      pcreadsend = __all_pos_creads[chr].second;

      npreadsstart = __all_neg_treads[chr].first;
      npreadsend = __all_neg_treads[chr].second;
      npcreadsstart = __all_neg_creads[chr].first;
      npcreadsend = __all_neg_creads[chr].second;
    } catch (...) {
      LOG_DEBUG1("CHR "<<chr<<" ignored. reads access exception.");
      break; //can not process this chr
    }
    //ignore the strand if it contains less than 2 read
    uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
    uint32_t nchrlength =
      npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));


    uint32_t noofbins = _calulate_number_of_bins(pchrlength,nchrlength,__binlength);
    uint32_t actualBinend = _calculate_last_pos(*(preadsend-1),*(npreadsend-1),__readlength,__readextension,pchrlength,nchrlength);
    uint32_t binind = 0;
    uint32_t binstart, binend;
    vector<called_peak> _result;
    _result.resize(0);
    vector<called_peak> regions;

    //TODO: Optimize so that it only process regions
    // that have reads

    LOG_DEBUG1("In Thread:" <<"Start processing " <<chr
        <<"\nNum of pos treat reads in this chr: "<<preadsend-preadsstart
        <<"\nThe first pos treat read in this chr: "<<*preadsstart
        <<"\nThe last pos treat read in this chr: "<<*(preadsend-1)
        <<"\npos length of this chr: "<<pchrlength
        <<"\nneg length of this chr: "<<nchrlength
        <<"\nNum of bins in this chr: "<<noofbins);

    //chromosomes -> pieces
    while (noofbins-- > 0) {

      /*
       * 1. generate a set of random reads and call peaks
       * 2. get the average num of peaks
       * 3. check if the current th gives
       *
       * fdr(threshold) = random peaks / actual peaks
       *
       * that fdr <= MIN_FDR
       *
       * 4. if 3 is true, return the actual peaks
       * otherwise, jump to 1 with th++.
       *
       *
       *
       */

      binstart = __binlength * binind + 1;
      binend = __binlength * (binind + 1);

      //peaks should not be larger than the last read
      if(binend >= actualBinend){
    	  binend = actualBinend;
    	  if(binend > 1){
    		  --binend; //in case of the 0-1 problem
    	  }
      }
      LOG_DEBUG1("binstart: "<<binstart<<"binend: "<<binend) ;
      reads_aux::reads_in(preadsstart, preadsend, binstart, binend,
          ppreadsstart, ppreadsend);
      LOG_DEBUG2("Num of pos treatment reads:"<<ppreadsend-ppreadsstart);
      reads_aux::reads_in(npreadsstart, npreadsend, binstart, binend,
          nppreadsstart, nppreadsend);
      LOG_DEBUG2("Num of neg treatment reads:"<<nppreadsend-nppreadsstart);

      if (ppreadsend == ppreadsstart && nppreadsend == nppreadsstart) {
        LOG_DEBUG2("NO reads are found for this bin. skip it");
        binind++;
        continue;
      }

      no_of_reads = (ppreadsend - ppreadsstart)
        + (nppreadsend - nppreadsstart);

      LOG_DEBUG1("Total reads:"<<no_of_reads);
      profile.resize(0);

      typedef vector<pair<uint32_t, int8_t> > reads_count_t;

      SGR sgr;
      vector<ActThresh> realPeaks;
      vector<Thresh> rndPeaks;
      uint32_t th = 0;

      // Generate SGR
      {
        get_profile_of_reads(__readextension, __readlength,
            ppreadsstart, ppreadsend, nppreadsstart, nppreadsend,
            sgr);
      } /* SGR generation */

      {

        findPeaksAtEachThreshold(sgr, realPeaks, binend);
        if (realPeaks.size() < 1)
          continue;
      } /* findpeaksatEachtreshold */

      if (realPeaks[0].peaks.size() < 1) {
        binind++;
        LOG_DEBUG2("NO real peaks were found for this bin. skip it");
        continue;
      }

      {
        getRndPeaksAtEachThreshold(binstart, binend, no_of_reads,
            __readextension, rndPeaks);

      } /* simulate */

      {
        th = findThresh(rndPeaks, realPeaks, fdr);
        LOG_DEBUG1("Estimated threshold: "<<th);
        assert(th >= MIN_THRESHOLD && th <= MAX_THRESHOLD);

      } /* Findthresh */

      regions.clear();
      //This is a must
      if (th >= MIN_THRESHOLD){
        th -= MIN_THRESHOLD;
      }
      else{
        cerr <<"bottom line threshold set as 1.\n";
        th = 1;
      }

      LOG_DEBUG2("realpeaks[th].peaks.size(): "<<realPeaks.at(th).peaks.size());
      assert_lt(th, realPeaks.size())
        regions = realPeaks.at(th).peaks;

      LOG_DEBUG2("Calculating p values. The configuration of thresholder:"
          << "\n__threshold:"<<th
          <<"\n__mergedistance:"<<200
          <<"\nNum of detected raw peaks:"<<regions.size());

      if (!(regions.empty())) {
        //if raw peaks were found, narrow the search range
        // of the control reads
        reads_aux::reads_in(pcreadsstart, pcreadsend, binstart, binend,
            ppcreadsstart, ppcreadsend);
        reads_aux::reads_in(npcreadsstart, npcreadsend, binstart,
            binend, nppcreadsstart, nppcreadsend);
      } else {
        LOG_DEBUG2("detected zero peaks for this bin. advance to the next bin ");
        binind++;
        continue;
      }

      //Get p value for each raw peak
      for (size_t i = 0; i < regions.size(); i++) {

        LOG_DEBUG3("Raw peak "<<i<<"\tpeak start: "
            <<regions[i].first<<"\tpeak end: "
            <<regions[i].second);
        if (regions[i].second > regions[i].first) {
          if (regions[i].second - regions[i].first < MIN_PEAKSIZE) {
            LOG_DEBUG3("Small peak, discarded");
            continue;
          }
          if (regions[i].second - regions[i].first > MAX_PEAKSIZE) {
            LOG_DEBUG3("Huge peak, discarded");
            continue;
          }
        } else {
          LOG_DEBUG3("peak end < peak start, discarded");
          continue;
        }
        //get reads from reads of the bin in each peak from both treatment and control
        reads_aux::reads_in(ppreadsstart, ppreadsend, regions[i].first,
            regions[i].second, pppreadsstart, pppreadsend);
        reads_aux::reads_in(ppcreadsstart, ppcreadsend,
            regions[i].first, regions[i].second, pppcreadsstart,
            pppcreadsend);
        reads_aux::reads_in(nppreadsstart, nppreadsend,
            regions[i].first, regions[i].second, npppreadsstart,
            npppreadsend);
        reads_aux::reads_in(nppcreadsstart, nppcreadsend,
            regions[i].first, regions[i].second, npppcreadsstart,
            npppcreadsend);

        a = (pppreadsend - pppreadsstart)
          + (npppreadsend - npppreadsstart);
        b = ceil(
            ((pppcreadsend - pppcreadsstart)
             + (npppcreadsend - npppcreadsstart))
            * norm_factor);

        LOG_DEBUG3("a "<<a);LOG_DEBUG3("b "<<b);

        try {
          p = distributions::binomCDF(a + b, b, 0.5);
        } catch (...) {
          printf(
              "Warning: possible overflow was caught in binom calculation."
              " a : %d\tb : %d\n", a, b);
          printf("pce-pcs:%d, npce-npcs:%d, norm_factor:%f\n",
              (uint32_t) (pppcreadsend - pppcreadsstart),
              (uint32_t) (npppcreadsend - npppcreadsstart),
              norm_factor);
          continue;
        }

        LOG_DEBUG3("p value: "<<p);

        if (p < __pvalThreshold) {

          if (print_stream) {
            vector<uint16_t> peak_region_profile;
            vector<uint32_t> max_pos, min_pos;

            try {
              try {
                region_profile::get_region_profile(
                    regions[i].first, regions[i].second,
                    __readextension, pppreadsstart,
                    pppreadsend, peak_region_profile);
              } catch (std::bad_alloc& e) {
               cerr 
                  << "Fatal: Not enough memory available while calulating "
                  "the profile of pos reads in summit detection. "
                  "The region in calculation is: ("
                  << regions[i].first << " : "
                  << regions[i].second << ")" << endl;

                exit(1);
              }
              try {
                region_profile::get_region_profile(
                    regions[i].first, regions[i].second,
                    __readextension, __readlength,
                    npppreadsstart, npppreadsend,
                    peak_region_profile, trans);
              } catch (std::bad_alloc& e) {
               cerr 
                  << "Fatal: Not enough memory available while calulating "
                  "the profile of neg reads in summit detection. "
                  "The region in calculation is: ("
                  << regions[i].first << " : "
                  << regions[i].second << endl;
                exit(1);
              }
              try {
                profile_smoother<uint16_t>::smooth(
                    peak_region_profile, __bandwidth);
              } catch (std::bad_alloc& e) {
               cerr 
                  << "Fatal: Not enough memory available while calulating "
                  "the smoothed profile in summit detection. "
                  "The region in calculation is: ("
                  << regions[i].first << " : "
                  << regions[i].second << endl;
                exit(1);
              }

              profile_summit_detector::detect_sub_peaks<uint16_t>(
                  peak_region_profile, __delta, max_pos,
                  min_pos);
            } catch (...) {
              LOG_DEBUG1( "Warning: Got errors while detecting summits" << endl);
              _profilenotgood = true;
            }
            boost::mutex::scoped_lock l(_result_mt);
            //todo: merge those boundary peaks here
            called_peak tmp;
            tmp.first = regions[i].first;
            tmp.second = regions[i].second;
            tmp.p = p;
            tmp.treads = a;
            tmp.creads = b;
            vector<uint32_t> tmpsumits;
            if (!_profilenotgood) {
              summits_2_vector<uint32_t>(max_pos,
                  regions[i].first, tmpsumits);
              tmp.summits = tmpsumits;
            }
            _result.push_back(tmp);
          }

          LOG_DEBUG3("Inserted peak "<<i
              <<"\tpeak start: "<<regions[i].first<<"\tpeak end: "
              <<regions[i].second);
        }
      }
      binind++; // advance to next bin
    }

    {
      boost::mutex::scoped_lock l(_result_mt);
      //consolidate discovered peaks and upload
      _resultRegions[chr] = _result;
      if (__verbose) {
        cerr << "Discovered " << _result.size() << "\tregions in "
          << chr << "." << endl;
      }
      // Fetech a new chr
      if (_workerPara._chrsleft.empty()) {
        //no more chrs to process
        break;//Only exit of the whole loop
      }
      chr = _workerPara._chrsleft.back();
      _workerPara._chrsleft.pop_back();
      // now using the internal _norm_factors
      norm_factor = _norm_factors[chr];
      LOG_DEBUG2("norm_factor:"<<chr<<":"<<norm_factor<<endl);
    }
  }
}

/*
 * Used by the wrapper
 */
void fdr_based_thresholder::detectSummits(Reads& treads, Reads& creads,
    uint32_t no_of_thread, double p_val_cutoff, double delta,
    uint32_t threshold, uint32_t mergedistance, uint32_t binlength,
    uint32_t bandwidth, uint32_t readextension, ostream& os, bool verbose) {

  //todo: treads.readlength == creads.readlength ?
  _setup(no_of_thread, binlength, threshold, mergedistance,
      treads.getReadlength(), readextension, bandwidth, p_val_cutoff,
      delta);

  _normalize_reads(treads, creads);
  /*
   * Divide the genome into chromosome
   * Call peaks in each chromosome by calling piceces of it
   * (now optional) consolidate peaks that span more than 1 regions
   * combine results
   */

  vector<uint32_t>::iterator treadsstart, treadsend, creadsstart, creadsend;
  /**
   * Each thread will fetch its own readsstart and readsend from this map
   */
  map<string, readsregion> all_pos_treads, all_neg_treads, all_pos_creads,
    all_neg_creads;

  /*
   * Prepare reads
   */
  assert_eq(treads.pos_reads.chrs().size(), treads.neg_reads.chrs().size());
  assert_eq(creads.pos_reads.chrs().size(), creads.neg_reads.chrs().size());
  assert_eq(treads.pos_reads.chrs().size(), creads.pos_reads.chrs().size());

  foreach(string chr, treads.pos_reads.chrs()) {
    treadsstart = treads.pos_reads.begin_of(chr);
    treadsend = treads.pos_reads.end_of(chr);
    creadsstart = creads.pos_reads.begin_of(chr);
    creadsend = creads.pos_reads.end_of(chr);
    all_pos_treads.insert(
        pair<string, readsregion>(chr,
          readsregion(treadsstart, treadsend)));
    all_pos_creads.insert(
        pair<string, readsregion>(chr,
          readsregion(creadsstart, creadsend)));
  }
  foreach(string chr, treads.neg_reads.chrs()) {
    treadsstart = treads.neg_reads.begin_of(chr);
    treadsend = treads.neg_reads.end_of(chr);
    creadsstart = creads.neg_reads.begin_of(chr);
    creadsend = creads.neg_reads.end_of(chr);
    all_neg_treads.insert(
        pair<string, readsregion>(chr,
          readsregion(treadsstart, treadsend)));
    all_neg_creads.insert(
        pair<string, readsregion>(chr,
          readsregion(creadsstart, creadsend)));
  }
  _workerPara.setAllPosTreads(all_pos_treads);
  _workerPara.setAllPosCreads(all_pos_creads);
  _workerPara.setAllNegTreads(all_neg_treads);
  _workerPara.setAllNegCreads(all_neg_creads);
  _workerPara.setVerbose(verbose);
  // Chrs are equal, assured by the assert and previous reads correction in peakranger.cpp
  _workerPara.setChrsleft(treads.pos_reads.chrs());

  //todo:* (now optional) consolidate peaks that span more than 1 regions

  boost::thread_group threads;
  size_t i = 0;
  while (i++ < _nThreads) {
    threads.create_thread(
        boost::bind(&fdr_based_thresholder::_processChr, this,
          boost::ref(os), true));

  }
  threads.join_all();
}

void fdr_based_thresholder::detectSummits(Reads& treads, Reads& creads,
    cmd_option_parser& option) {

  ofstream os;
  detectSummits(treads, creads, option.getNo_of_thread(), option.getCut_off(),
      option.getDelta(), option.getPeakHeightCutoff(), MERGE_DISTANCE,
      option.getBinlength(), option.getBandwidth(),
      option.getExt_length(), os, option.getVerboseRequested());
}
/*
 * Wrapper, the entry function used
 */
void fdr_based_thresholder::detectSummits(Reads & treads, Reads & creads,
    cmd_option_parser & option, ostream & os) {
  os << "#  Results generated by PeakRanger.\n"
    "#  Please cite:\n"
    "#  Feng X, Grossman R, Stein L: PeakRanger: \n"
    "#  A cloud-enabled peak caller for ChIP-seq data. \n"
    "#  BMC Bioinformatics 2011, 12(1):139.\n\n#";
  os
    << "#chr\tstart\tend\tsummits\tvalleys\tp-value\tsample_reads\tcontrol_reads\n";
  LOG_DEBUG("\nParameters for simple thresholder detector:"
      <<"\nthreads:"<< option.getNo_of_thread()
      <<"\np value cut off:"<<option.getCut_off()
      <<"\ndelta:"<<option.getDelta()
      <<"\npeak height cutoff:"<<option.getPeakHeightCutoff()
      <<"\nmerge distance:"<<MERGE_DISTANCE
      <<"\nbin length:"<<option.getBinlength()
      <<"\nbandwidth:"<<option.getBandwidth()
      <<"\nextension length:"<<option.getExt_length());

  detectSummits(treads, creads, option.getNo_of_thread(), option.getCut_off(),
      option.getDelta(), option.getPeakHeightCutoff(), MERGE_DISTANCE,
      option.getBinlength(), option.getBandwidth(),
      option.getExt_length(), os, option.getVerboseRequested());

}

/*
 * Wrapper for normalization.
 *
 * treads and creads must contain at least one chr
 */
void fdr_based_thresholder::_normalize_reads(Reads& treads, Reads& creads) {

  assert_eq(treads.pos_reads.chrs().size(), treads.neg_reads.chrs().size());
  assert_eq(creads.pos_reads.chrs().size(), creads.neg_reads.chrs().size());
  assert_eq(treads.pos_reads.chrs().size(), creads.pos_reads.chrs().size());

  foreach(string chr, treads.pos_reads.chrs()) {
    _norm_factors[chr] = _normalize_reads_onchr(treads, creads, chr,
        _binlength);

  }
  foreach(string chr, treads.pos_reads.chrs()) {
    if (_norm_factors[chr] == 1) {
      LOG_DEBUG2( "Warning: Normalization for " << chr << " failed." << endl);
    }
  }
}

/*  Calculate the normalization factor of creads toward treads.
 *
 * todo: treads toward creads?
 *
 * pre-requisite:
 * return: a value ~0 ~ 1. If the normalization cant finish, return 1
 */
inline double fdr_based_thresholder::_normalize_reads_onchr(Reads& treads,
    Reads& creads, string& chr, uint32_t binlength) {
  /*
   * 	Least square regression
   *
   *   (Sigma(xy) - n*MEAN(x)MEAN(y))
   * ----------------------------------
   *   (Signma(x^2) - n*((MEAN(x))^2))
   *
   *   http://mathworld.wolfram.com/LeastSquaresFitting.html
   */

  double result = 1;

  vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
    preadsend, ppreadsend, pppreadsend, pcreadsstart, ppcreadsstart,
    pppcreadsstart, pcreadsend, ppcreadsend, pppcreadsend;

  vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
    npreadsend, nppreadsend, npppreadsend, npcreadsstart,
    nppcreadsstart, npppcreadsstart, npcreadsend, nppcreadsend,
    npppcreadsend;

  try {
    preadsstart = treads.pos_reads.begin_of(chr);
    preadsend = treads.pos_reads.end_of(chr);
    npreadsstart = treads.neg_reads.begin_of(chr);
    npreadsend = treads.neg_reads.end_of(chr);
    pcreadsstart = creads.pos_reads.begin_of(chr);
    pcreadsend = creads.pos_reads.end_of(chr);
    npcreadsstart = creads.neg_reads.begin_of(chr);
    npcreadsend = creads.neg_reads.end_of(chr);
  } catch (...) {
    return 1;
  }

  uint64_t pchrlength = (*(preadsend - 1));
  uint64_t nchrlength = (*(npreadsend - 1));
  uint64_t apchrlength = pchrlength > nchrlength ? pchrlength : nchrlength;
  /*
   *
   * todo: Must verify the bin length now to avoid overflow
   *
   * binlength^2 <= 0xffffffff
   */
  if (binlength > 100000) {
    cerr << "binlength is larger than " << 100000 << endl;
    binlength = 0x3fff;
    cerr << "changed to " << binlength << endl;
  }
  uint64_t noofbins = apchrlength ? 1 + (apchrlength / binlength) : 0;
  uint64_t binind = 0, actualbins = 0;
  uint64_t binstart, binend;
  uint64_t treadsno = 0, creadsno = 0, streadsno = 0, screadsno = 0,
           sccreadsno = 0, stcreadsno = 0;
  uint64_t t_pre_pos = 0, nt_pre_pos = 0, c_pre_pos = 0, nc_pre_pos = 0;
  //todo: how to effectively remove reads in peak regions?
  //response: the following codes generated satisfactory factors
  //todo: reads with the same start positions should be counted only once
  //response: done, see the codes.

  LOG_DEBUG1("\n\n*******Start calculating scaling factor for "<<chr<<"******"
      <<"\nNo of bins : "<<noofbins<<"\tBin size : "<<binlength);
  while (noofbins-- > 0) {
    treadsno = 0, creadsno = 0;
    binstart = binlength * binind + 1;
    binend = binlength * (binind + 1);
    LOG_DEBUG2("Processing bin : "<<binind<<
        "\nbin_start : "<<binstart<<"\tbin_stop : "<<binend);
    ppreadsstart = lower_bound(preadsstart, preadsend, binstart);
    ppreadsend = upper_bound(ppreadsstart, preadsend, binend);
    ppcreadsstart = lower_bound(pcreadsstart, pcreadsend, binstart);
    ppcreadsend = upper_bound(ppcreadsstart, pcreadsend, binend);

    nppreadsstart = lower_bound(npreadsstart, npreadsend, binstart);
    nppreadsend = upper_bound(nppreadsstart, npreadsend, binend);
    nppcreadsstart = lower_bound(npcreadsstart, npcreadsend, binstart);
    nppcreadsend = upper_bound(nppcreadsstart, npcreadsend, binend);

    while (ppcreadsstart != ppcreadsend) {
      if (*ppcreadsstart == c_pre_pos) {
        ppcreadsstart++;
        continue;
      }
      creadsno++;
      c_pre_pos = *ppcreadsstart;
      ppcreadsstart++;
    }
    while (nppcreadsstart != nppcreadsend) {
      if (*nppcreadsstart == nc_pre_pos) {
        nppcreadsstart++;
        continue;
      }
      creadsno++;
      nc_pre_pos = *nppcreadsstart;
      nppcreadsstart++;
    }
    if (creadsno == 0) {
      binind++;
      continue;
    }
    while (ppreadsstart != ppreadsend) {
      if (*ppreadsstart == t_pre_pos) {
        ppreadsstart++;
        continue;
      }
      treadsno++;
      t_pre_pos = *ppreadsstart;
      ppreadsstart++;
    }
    while (nppreadsstart != nppreadsend) {
      if (*nppreadsstart == nt_pre_pos) {
        nppreadsstart++;
        continue;
      }
      treadsno++;
      nt_pre_pos = *nppreadsstart;
      nppreadsstart++;
    }
    if (treadsno == 0) {
      binind++;
      continue;
    }
    streadsno += treadsno;
    screadsno += creadsno;
    // todo: overflow?
    stcreadsno += treadsno * creadsno;
    sccreadsno += creadsno * creadsno;
    LOG_DEBUG2("total control reads of this bin:"<<creadsno);LOG_DEBUG2("total sample reads of this bin:"<<treadsno);
    binind++;
    actualbins++;
  }LOG_DEBUG1("final values for this bin:");LOG_DEBUG1("n : " << actualbins);LOG_DEBUG1("sxy : " << stcreadsno);LOG_DEBUG1("sxx : " << sccreadsno);LOG_DEBUG1("sx : " << screadsno);LOG_DEBUG1("sy : " << streadsno);

  if (actualbins == 1) {
    return 1;
  }

  if (actualbins * stcreadsno > streadsno * screadsno) {

    if (actualbins * sccreadsno > screadsno * screadsno) {
      result = (actualbins * stcreadsno - streadsno * screadsno) * 1.0
        / (actualbins * sccreadsno - screadsno * screadsno);
      if (result < 0) {
        return 1;
      }LOG_DEBUG1("final value : " << (actualbins * stcreadsno - streadsno * screadsno) * 1.0
          / (actualbins * sccreadsno - screadsno * screadsno));
    }
  }

  return result;

}

map<string, readsregion> fdr_based_thresholder::workerPara::getAllNegCreads() const {
  return _all_neg_creads;
}

map<string, readsregion> fdr_based_thresholder::workerPara::getAllNegTreads() const {
  return _all_neg_treads;
}

map<string, readsregion> fdr_based_thresholder::workerPara::getAllPosCreads() const {
  return _all_pos_creads;
}

map<string, readsregion> fdr_based_thresholder::workerPara::getAllPosTreads() const {
  return _all_pos_treads;
}

vector<string> fdr_based_thresholder::workerPara::getChrsleft() const {
  return _chrsleft;
}

bool fdr_based_thresholder::workerPara::isVerbose() const {
  return _verbose;
}

void fdr_based_thresholder::workerPara::setAllNegCreads(
    map<string, readsregion> _all_neg_creads) {
  this->_all_neg_creads = _all_neg_creads;
}

void fdr_based_thresholder::workerPara::setAllNegTreads(
    map<string, readsregion> _all_neg_treads) {
  this->_all_neg_treads = _all_neg_treads;
}

void fdr_based_thresholder::workerPara::setAllPosCreads(
    map<string, readsregion> _all_pos_creads) {
  this->_all_pos_creads = _all_pos_creads;
}

void fdr_based_thresholder::workerPara::setAllPosTreads(
    map<string, readsregion> _all_pos_treads) {
  this->_all_pos_treads = _all_pos_treads;
}

void fdr_based_thresholder::workerPara::setChrsleft(vector<string> _chrsleft) {
  this->_chrsleft = _chrsleft;
}

void fdr_based_thresholder::workerPara::setVerbose(bool _verbose) {
  this->_verbose = _verbose;
}

uint32_t fdr_based_thresholder::_calulate_number_of_bins(uint32_t p_chr_length, uint32_t n_chr_length,
		uint32_t binlength) const{

	 uint32_t apchrlength =
			 p_chr_length > n_chr_length ? p_chr_length : n_chr_length;
	 if (apchrlength > 0) {
		 if (apchrlength < binlength) {
			 return 1;
		 }
		 if (apchrlength == binlength){
			 return  1;
		 }
		 if (apchrlength > binlength){
			 if(apchrlength % binlength == 0){
				 return apchrlength / binlength;
			 }
			 else{
				 return 1+(apchrlength / binlength);
			 }
		 }
	 }

		 return 0;

}

uint32_t fdr_based_thresholder::_calculate_neg_read_loci_ret_start(
		uint32_t start, uint32_t readlength, uint32_t extension) const{
	 return start + readlength - extension;
}

uint32_t fdr_based_thresholder::_calculate_neg_read_loci_ret_end(uint32_t start,
		uint32_t readlength, uint32_t extension) const{
	return start + readlength; //changed to match JTwiggle
}

uint32_t fdr_based_thresholder::_calculate_pos_read_loci_ret_start(
		uint32_t start, uint32_t readlength, uint32_t extension) const{
	return start;
}

uint32_t fdr_based_thresholder::_calculate_pos_read_loci_ret_end(uint32_t start,
		uint32_t readlength, uint32_t extension) const{

	return start+extension;
}
