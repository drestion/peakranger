/*
 * ccat_aux.cpp
 *
 *  Created on: Mar 28, 2012
 *      Author: xin
 */
#include "ccat_aux.h"
#include "common/stl_header.h"
#include "common/ranger_debug.h"
#include "common/boost_header.h"
#include "short_reads/readstools.h"
using namespace std;
namespace ccat_aux {
  void load_histone_config(ccat_config_t & config) {
    config.fragmentSize = 200;
    config.slidingWinSize = 500;
    config.movingStep = 50;
    config.isStrandSensitiveMode = 0;
    config.minCount = 4;
    config.outputNum = 100000;
    config.randomSeed = 123456;
    config.minScore = 5.0;
    config.bootstrapPass = 50;
  }

  void config_validation(ccat_config_t& config) {
    if (config.fragmentSize < 0) {
      config.fragmentSize = 200;
    }
    if (config.slidingWinSize < 1) {
      config.slidingWinSize = 300;
    }
    if (config.movingStep < 1) {
      config.movingStep = 10;
    }
    if (config.outputNum < 0) {
      config.outputNum = 100000;
    }
    if (config.minCount <= 1) {
      config.minCount = 2;
    }
    if (config.minScore < 0) {
      config.minScore = 5.0;
    }
    if (config.bootstrapPass < 1) {
      config.bootstrapPass = 50;
    }
    if (config.bootstrapPass > 50) {
      config.bootstrapPass = 50;
    }
  }

  void GenRegionProfile(vector<size_t>& l1Profile, vector<size_t>& l2Profile,
      vector<size_t>& region, double l1Ratio, double l2Ratio,
      size_t slidingWinSize, size_t movingStep, size_t bootstrapPass,
      double minScore, size_t row, size_t column, size_t tagCount,
      size_t binCount, double smoothingFactor, vector<double>& lookUpTable,
      vector<int>& flag) {

    int i, j;
    double tmpF;
    int halfWinSize = slidingWinSize / movingStep / 2;
    size_t len = l1Profile.size();
    assert(len > 1);
    assert(l1Profile.size() == l2Profile.size());
    assert(l1Profile.size() == region.size());
    vector<size_t> tmpRegion(len, 0);
    for (i = 0; i < len; i++) {
      if ((l1Profile[i] >= row) || (l2Profile[i] >= column)) {

        throw RangerException(
            "In GenRegionProfile: l1Profile or l2Profile is larger than the max possible value.\n");
      }
      if (flag[l1Profile[i] * column + l2Profile[i]]) {
        tmpRegion[i] =
          lookUpTable[l1Profile[i] * column + l2Profile[i]]
          > minScore ? 1 : 0;
        continue;
      }
      tmpF = BootstrapFoldChange(l1Profile[i], l2Profile[i], l1Ratio, l2Ratio,
          bootstrapPass, tagCount, binCount, smoothingFactor);

      if (tmpF > minScore) {
        tmpRegion[i] = 1;
      } else {
        tmpRegion[i] = 0;
      }
      lookUpTable[l1Profile[i] * column + l2Profile[i]] = tmpF;
      flag[l1Profile[i] * column + l2Profile[i]] = 1;
    }

    assert(region.size() >= tmpRegion.size());
    copy(tmpRegion.begin(), tmpRegion.end(), region.begin());
    for (i = 1; i < len - 1; i++) {
      if ((tmpRegion[i] == 1) && (tmpRegion[i - 1] == 0)) {
        for (j = i - 1; j >= i - halfWinSize; j--) {
          if ((tmpRegion[j]) || (j < 0)) {
            break;
          }
          region[j] = 1;
        }

      }

      if ((tmpRegion[i] == 1) && (tmpRegion[i + 1] == 0)) {
        for (j = i + 1; j <= i + halfWinSize; j++) {
          if ((tmpRegion[j]) || (j >= len)) {
            break;
          }
          region[j] = 1;
        }

      }

    }

  }

  double BootstrapFoldChange(size_t l1Count, size_t l2Count, double l1Ratio,
      double l2Ratio, size_t bootstrapPass, size_t tagCount, size_t binCount,
      double smoothingFactor) {
    int count1, count2;
    double foldChange;
    vector<double> tmpA(bootstrapPass, 0);
    size_t i, j;
    int ga;

    for (i = 0; i < bootstrapPass; i++) {
      count1 = 0;

      for (j = 0; j < l1Count; j++) {
        ga = rand();
        if (ga > RAND_MAX * l1Ratio) {
          continue;
        }

        count1++;
      }

      count2 = 0;

      for (j = 0; j < l2Count; j++) {
        ga = rand();
        if (ga > RAND_MAX * l2Ratio) {
          continue;
        }

        count2++;
      }

      if (count1 + count2 == 0) {
        tmpA[i] = 0;
      } else {

        tmpA[i] = ComputeFoldChange(count1, count2, tagCount, binCount,
            smoothingFactor);

      }
    }
    foldChange = 0.0;
    for (i = 0; i < bootstrapPass; i++) {
      foldChange += tmpA[i];
    }
    return foldChange / bootstrapPass;
  }

  double ComputeFoldChange(size_t count1, size_t count2, size_t tagCount,
      size_t binCount, double smoothingFactor) {
    return (double) (count1) / (count2 + smoothingFactor)
      * (tagCount + smoothingFactor * binCount) / tagCount;
  }

  double ComputeSmoothingParameter(vector<peak_t>& peaks, size_t& binCount,
      size_t& tagCount) {
    size_t hist_size = 1000;
    size_t peakNum = peaks.size();
    vector<int> hist(hist_size, 0);
    double m=0, result=0;
    for (size_t i = 0; i < peakNum; i++) {
      if (peaks[i].reSampledL2Count < 1000) {
        hist[peaks[i].reSampledL2Count]++;

        binCount++;
        tagCount += peaks[i].reSampledL2Count;
      }
    }

    FitNegBinomDist(hist, m, result, 0.1, hist_size * 1.0);
    return result;
  }

  int FitNegBinomDist(vector<int>& hist, double& m, double& k, double k_min,
      double k_max) {
    size_t i, j;
    double sum;
    double tmp_k_min, tmp_k_max, tmpK, step;
    double error, minError;
    double bestFit=0, prevBestFit=0;
    size_t histLen = hist.size();
    vector<double> histNorm(histLen, 0);
    if ((k_max <= k_min) || (k_min <= 0)) {
      return -1;
    }
    sum = 0.0;
    foreach(int h, hist) {
      sum += h;
    }
    if (sum < 1.0) {
      return -1;
    }
    m = 0.0;
    for (i = 0; i < histLen; i++) {
      histNorm[i] = (double) ((hist[i])) / sum;
      m += histNorm[i] * (double) (i);
    }
    tmp_k_min = k_min;
    tmp_k_max = k_max;
    prevBestFit = DBL_MAX;
    while (1) {
      minError = DBL_MAX;

      step = log(tmp_k_max / tmp_k_min) / 9;

      for (i = 0; i < 10; i++) {
        tmpK = tmp_k_min * exp(step * (double) i);

        error = 0.0;

        sum = 0.0;

        for (j = 0; j < histLen; j++) {
          error += (exp(
                gammln(tmpK + j) - gammln(tmpK) - gammln((double) j + 1)
                + tmpK * log(tmpK / (tmpK + m))
                + (double) j * log(m / (tmpK + m)))
              - histNorm[j])
            * (exp(
                  gammln(tmpK + j) - gammln(tmpK)
                  - gammln((double) j + 1)
                  + tmpK * log(tmpK / (tmpK + m))
                  + (double) j * log(m / (tmpK + m)))
                - histNorm[j]);
        }
        error = sqrt(error);

        if (error < minError) {
          bestFit = tmpK;
          minError = error;
        }
      }

      if (fabs(log(prevBestFit / bestFit)) < 0.01) {
        k = bestFit;
        break;
      }

      tmp_k_min = bestFit / exp(step) < k_min ? k_min : bestFit / exp(step);
      tmp_k_max = bestFit * exp(step) > k_max ? k_max : bestFit * exp(step);

      prevBestFit = bestFit;
    }
    return 1;
  }

  void initializeSrandUsingCurrentTime() {
    srand(time(NULL));
  }

  double gammln(double xx) {
    double x, y, tmp, ser;
    double cof[6] = { };
    int j;

    if (xx < 0.0) {
      return -1.0;
    }

    if (xx == 0.0) {
      return 0.0;
    }

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;

    for (j = 0; j <= 5; j++) {
      ser += cof[j] / ++y;
    }

    return -tmp + log(2.5066282746310005 * ser / x);
  }

  void LoadData(vector<chr_t>& chroms, Reads& treads, Reads& creads,
      size_t& chromNum) {

    assert_eq(treads.pos_reads.chrs().size(), treads.neg_reads.chrs().size());
    assert_eq(creads.pos_reads.chrs().size(), creads.neg_reads.chrs().size());
    assert_eq(treads.pos_reads.chrs().size(), creads.pos_reads.chrs().size());
    chromNum = treads.pos_reads.chrs().size();

    size_t i = 0;
    foreach(string chr, treads.pos_reads.chrs()) {
      //todo: assuming reads are sorted.
      chr_t chrom;
      chrom.chromName = chr;
      chrom.chromIndex = i;
      size_t treads_chr_size = reads_tools::chromSize(treads, chr);
      size_t creads_chr_size = reads_tools::chromSize(creads, chr);

      chrom.chromSize =
        treads_chr_size > creads_chr_size ?
        treads_chr_size : creads_chr_size;

      copy(treads.pos_reads.begin_of(chr), treads.pos_reads.end_of(chr),
          std::back_inserter(chrom.l1PosTags));
      copy(treads.neg_reads.begin_of(chr), treads.neg_reads.end_of(chr),
          std::back_inserter(chrom.l1NegTags));
      copy(creads.pos_reads.begin_of(chr), creads.pos_reads.end_of(chr),
          std::back_inserter(chrom.l2PosTags));
      copy(creads.neg_reads.begin_of(chr), creads.neg_reads.end_of(chr),
          std::back_inserter(chrom.l2NegTags));

      chroms.push_back(chrom);
      i++;
    }

  }

  void SortAndDedup(vector<chr_t>&chroms) {

    for (size_t i = 0; i < chroms.size(); i++) {
      ccat_aux::SortAndGetUniqueTags(chroms[i].l1PosTags);
      ccat_aux::SortAndGetUniqueTags(chroms[i].l2PosTags);
      ccat_aux::SortAndGetUniqueTags(chroms[i].l1NegTags);
      ccat_aux::SortAndGetUniqueTags(chroms[i].l2NegTags);
    }
  }
} //namespace ccat_aux
