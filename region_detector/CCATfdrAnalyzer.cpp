/*
 * CCATfdrAnalyzer.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#include "CCATfdrAnalyzer.h"
#include "common/ranger_debug.h"
#include "ccat_aux.h"
#include "common/boost_header.h"
#include "utils/stl_helper.h"
#include "ccat_profile_helper.h"
using namespace std;
namespace ccat_aux {
  const size_t CCATfdrAnalyzer::QVALUESTEP = 1000;
  CCATfdrAnalyzer::CCATfdrAnalyzer() :
    q(QVALUESTEP, 0), value(QVALUESTEP, 0), tagCount(0), binCount(0) {

    }

  CCATfdrAnalyzer::~CCATfdrAnalyzer() {

  }

  //SignificanceAnalysis: perform significance analysis, will be called by the main routine of CCAT
  int CCATfdrAnalyzer::SignificanceAnalysis(vector<chr_t>&chroms, size_t chromNum, double l1Ratio, double l2Ratio,
      size_t maxL1Count, size_t maxL2Count, const ccat_config_t & config) {
    MARK_FUN("CCATfdrAnalyzer::SignificanceAnalysis");
    size_t l1PeakNum, l2PeakNum;
    size_t i, j;
    double threshFDR;

    l1PeakNum = 0;
    l2PeakNum = 0;

    for (i = 0; i < chromNum; i++) {
      l1PeakNum += chroms[i].l1Peaks.size();
      l2PeakNum += chroms[i].l2Peaks.size();
    }

    if ((l1PeakNum <= 0) || (l2PeakNum <= 0)) {
      return -1;
    }

    vector<peak_t> l1Peaks(l1PeakNum);
    vector<peak_t> l2Peaks(l2PeakNum);
    l1PeakNum = 0;
    l2PeakNum = 0;

    //re-sampling

    for (i = 0; i < chromNum; i++) {
      for (j = 0; j < chroms[i].l1Peaks.size(); j++) {

        l1Peaks[l1PeakNum] = chroms[i].l1Peaks[j];
        l1PeakNum++;
      }

      for (j = 0; j < chroms[i].l2Peaks.size(); j++) {
        l2Peaks[l2PeakNum] = chroms[i].l2Peaks[j];
        l2PeakNum++;
      }
    }
    smoothingFactor = ComputeSmoothingParameter(l1Peaks, binCount, tagCount);
    if (smoothingFactor < 0) {
      smoothingFactor = 1.0;
    }
    threshFDR = ComputeLocalFDR(l1Peaks, l2Peaks, q, value, config, tagCount, binCount);
    if (threshFDR > 0.2) {
      //todo: this must not be set if I want my opt results.
      //        cout << "Warning: estimated cutoff FDR is "<<threshFDR<<".\n");
    } else {
      cout << ("\n");
    }

    PostProcessing(chroms, chromNum, l1Ratio, l2Ratio, maxL1Count, maxL2Count, config);
    return 1;
  }

  //ComputeThreshold: compute the local FDR, return the FDR at the cut-off threshold;
  double CCATfdrAnalyzer::ComputeLocalFDR(vector<peak_t>& l1Peaks, vector<peak_t>& l2Peaks, vector<double>&q,
      vector<double> &value, const ccat_config_t& config, size_t tagCount, size_t binCount) {
    MARK_FUN("CCATfdrAnalyzer::ComputeLocalFDR");

    LOG_DEBUG2("tagCount:"<<tagCount);LOG_DEBUG2("binCount:"<<binCount);

    LOG_DEBUG2("smoothingFactor:"<<smoothingFactor);

    size_t l1PeakNum = l1Peaks.size();
    size_t l2PeakNum = l2Peaks.size();

    LOG_DEBUG2("l1PeakNum:"<<l1PeakNum);LOG_DEBUG2("l2PeakNum:"<<l2PeakNum);

    size_t posCount, negCount;
    size_t tmpIndex;
    int i; //cant be size_t
    size_t j;
    double score;

    vector<double> p1(l1PeakNum, 0);
    vector<double> p2(l2PeakNum, 0);
    double threshFDR = 0;

    if (l1PeakNum < 1 || l2PeakNum < 1) {
      return threshFDR;
    }

    for (i = 0; i < l1PeakNum; i++) {
      p1[i] = ComputeFoldChange(l1Peaks[i].reSampledL1Count, l1Peaks[i].reSampledL2Count, tagCount, binCount,
          smoothingFactor);
    }

    for (i = 0; i < l2PeakNum; i++) {
      p2[i] = ComputeFoldChange(l2Peaks[i].reSampledL1Count, l2Peaks[i].reSampledL2Count, tagCount, binCount,
          smoothingFactor);
    }

    std::sort(p1.begin(), p1.end(), greater<double>());
    std::sort(p2.begin(), p2.end(), greater<double>());

    for (j = 0; j < QVALUESTEP; j++) {
      q[j] = (double) (j + 1) / QVALUESTEP;

      score = 0;
      tmpIndex = -1;

      for (i = static_cast<int>(l1PeakNum - 1); i >= 0; i--) {
        if (i < l1PeakNum - 1) {
          if (p1[i] - p1[i + 1] < 0.00000001) {
            continue;
          }
        }
        posCount = i + 1;
        negCount = lower_bound(p2.begin(), p2.end(), p1.at(i) - 0.00000001, greater<double>()) - p2.begin();
        if (posCount * q[j] - negCount > score) {
          score = posCount * q[j] - negCount;
          tmpIndex = i;
        }
      }

      if (tmpIndex == -1) {
        value[j] = DBL_MAX;
      } else {
        value[j] = p1[tmpIndex];
      }
    }

    threshFDR = 0.0;

    for (i = 1; i < QVALUESTEP; i++) {
      if (value[i] > value[i - 1]) {
        value[i] = value[i - 1];
      }

      if (value[i] > config.minScore) {
        threshFDR = (double) i / QVALUESTEP;
      }
    }
    return threshFDR;
  }

  //PostProcessing: post-processing of the peaks: boostrapping for fold-change calculation, region identification, peak refinement
  int CCATfdrAnalyzer::PostProcessing(vector<chr_t>& chroms, int chromNum, double l1Ratio, double l2Ratio, int maxL1Count,
      int maxL2Count, const ccat_config_t& config) {
    MARK_FUN("CCATfdrAnalyzer::PostProcessing");
    int i;
    row = maxL1Count + 1;
    column = maxL2Count + 1;
    LOG_DEBUG2("Row "<< row);LOG_DEBUG2("Column "<< column);
    size_t size = row * column;
    vector<double> tmp(size, 0);
    vector<int> tp(size, 0);
    lookUpTable.swap(tmp);
    flag.swap(tp);

    for (i = 0; i < chromNum; i++) {
      if (ProcessOneChrom(chroms[i], l1Ratio, l2Ratio, config) < 0) {
        return -1;
      }
    }

    return 1;
  }

  //ProcessOneChrom: process one chromosome
  int CCATfdrAnalyzer::ProcessOneChrom(chr_t& chrom, double l1Ratio, double l2Ratio, const ccat_config_t& config) {
    MARK_FUN("CCATfdrAnalyzer::ProcessOneChrom");
    size_t i;
    int profileLen = (chrom.chromSize) / config.movingStep + 1;
    int tmpIndex, tmpStart, tmpEnd, lastStart, lastEnd;
    vector<size_t> region(profileLen, 0);
    vector<size_t> profile1(region);
    vector<size_t> profile2(region);
    //generate profile, assign reads from both ChIP and control libraries to the profile.
    profile.setFragmentSize(config.fragmentSize);
    profile.setMovingStep(config.movingStep);
    profile.setProfileLength(profileLen);
    buildL1Profile(profile, chrom, profile1);
    buildL2Profile(profile, chrom, profile2);
    profile.SmoothProfile(profile1, profileLen, config.slidingWinSize / config.movingStep / 2);
    profile.SmoothProfile(profile2, profileLen, config.slidingWinSize / config.movingStep / 2);

    GenRegionProfile(profile1, profile2, region, l1Ratio, l2Ratio, config.slidingWinSize, config.movingStep,
        config.bootstrapPass, config.minScore, row, column, tagCount, binCount, smoothingFactor, lookUpTable, flag);

    lastStart = -1;
    lastEnd = -1;
    for (i = 0; i < chrom.l1Peaks.size(); i++) {
      chrom.l1Peaks[i].foldChange = lookUpTable[chrom.l1Peaks[i].l1Count * column + chrom.l1Peaks[i].l2Count];
      LOG_DEBUG4(" chrom.l1Peaks[i].foldChange = "<< chrom.l1Peaks[i].foldChange);

      tmpIndex = lower_bound(value.begin(), value.end(), chrom.l1Peaks[i].foldChange - 0.00000001, greater<double>())
        - value.begin();
      if (tmpIndex != 0) {
        tmpIndex -= 1; // -1 is matched for the original bTreeSort.
      }

      //Measured significance is foldChange
      chrom.l1Peaks[i].qValue = q[tmpIndex];
      LOG_DEBUG4("chrom.l1Peaks["<<i<<"].qValue = q["<<tmpIndex<<"];");
      if (chrom.l1Peaks[i].foldChange > config.minScore) {
        chrom.l1Peaks[i].isSignificant = 1;
      } else {
        chrom.l1Peaks[i].isSignificant = 0;
      }

      if (chrom.l1Peaks[i].isSignificant) {
        if (chrom.l1Peaks[i].peak <= lastEnd) {
          tmpStart = lastStart;
          tmpEnd = lastEnd;
        } else {
          for (tmpStart = chrom.l1Peaks[i].peak - 1; tmpStart >= 0; tmpStart--) {
            if (!region[tmpStart]) {
              break;
            }
          }

          tmpStart++;

          for (tmpEnd = chrom.l1Peaks[i].peak + 1; tmpEnd < profileLen; tmpEnd++) {
            if (!region[tmpEnd]) {
              break;
            }
          }

          tmpEnd--;

          lastStart = tmpStart;
          lastEnd = tmpEnd;
        }

        chrom.l1Peaks[i].start = tmpStart * config.movingStep;
        chrom.l1Peaks[i].end = (tmpEnd + 1) * config.movingStep;
      } else {
        chrom.l1Peaks[i].start = chrom.l1Peaks[i].peak * config.movingStep + config.movingStep / 2;
        chrom.l1Peaks[i].end = chrom.l1Peaks[i].peak * config.movingStep + config.movingStep / 2;
      }

      chrom.l1Peaks[i].peak = chrom.l1Peaks[i].peak * config.movingStep + config.movingStep / 2;
    }
    return 1;
  }

} /* namespace ccat_aux */
