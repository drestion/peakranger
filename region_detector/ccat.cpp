/*
 * ccat.cpp
 *
 *  Created on: Dec 14, 2011
 *      Author: xfeng
 */
#include "ccat.h"
#include "common/boost_header.h"
#include "common/stl_header.h"
#include "utils/assert_helpers.h"
#include "utils/logger.h"
#include "utils/debug.h"
#include "short_reads/readstools.h"
#include "region_detector/calledpeak.h"

using namespace std;
using namespace ccat_aux;

ccat::ccat() :
  lookUpTable(), flag(), q(Q_VALUE_STEP_NUM, 0), value(Q_VALUE_STEP_NUM,
      0), conf(), prof(), noise(), fdr(), pkFinder() {

  }

ccat::~ccat() {

}

void ccat::LoadData(vector<chr_t>& chroms, Reads& treads, Reads& creads,
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
void ccat::insertPeak(const string& chr, called_peak& pk) {
  map<string, vector<called_peak> >::iterator it;
  it = _resultRegions.find(chr);

  if (it == _resultRegions.end()) {
    vector<called_peak> peaks;
    _resultRegions[chr] = peaks;
  }
  _resultRegions[chr].push_back(pk);
}

//OutputFiles: output the results
int ccat::UploadPeaks(vector<chr_t>&chroms, size_t chromNum,
    const char *projectName, const ccat_config_t& config) {
  size_t peakNum;
  size_t i, j;

  string fileName(projectName);
  size_t startChrom;
  int sumL1, sumL2;

  peakNum = 0;
  sumL1 = 0;
  sumL2 = 0;

  for (i = 0; i < chromNum; i++) {
    peakNum += chroms[i].l1Peaks.size();
    sumL1 += chroms[i].l1PosTags.size() + chroms[i].l1NegTags.size();
    sumL2 += chroms[i].l2PosTags.size() + chroms[i].l2NegTags.size();
  }

  if (peakNum <= 0) {
    return 0;
  }

  vector<peak_t> peaks(peakNum);

  peakNum = 0;

  for (i = 0; i < chromNum; i++) {

    for (j = 0; j < chroms[i].l1Peaks.size(); j++) {

      peaks[peakNum] = chroms[i].l1Peaks[j];
      peaks[peakNum].chromIndex = i;
      peakNum++;
    }
  }

  sort(peaks.begin(), peaks.end(), greaterFoldChange);

  for (startChrom = 0; startChrom < chromNum; startChrom++) {
    if (chroms[startChrom].l1Peaks.size() > 0) {
      break;
    }
  }

  peaks[0] = chroms[startChrom].l1Peaks[0];
  peaks[0].chromIndex = startChrom;
  peakNum = 1;
  for (i = startChrom; i < chromNum; i++) {
#ifdef DEBUG
    string ff("fc");
    ff += chroms[i].chromName;
    ranger_debug::dumpArray<vector<peak_t>, ccat_aux::mm>(peaks, ff.c_str());
#endif
    for (j = 0; j < chroms[i].l1Peaks.size(); j++) {
      if ((i == startChrom) && (j == 0)) {
        continue;
      }

      if ((i == peaks[peakNum - 1].chromIndex)
          && (chroms[i].l1Peaks[j].start == peaks[peakNum - 1].start)) {
        if (chroms[i].l1Peaks[j].foldChange
            > peaks[peakNum - 1].foldChange) {
          peaks[peakNum - 1] = chroms[i].l1Peaks[j];
          peaks[peakNum - 1].chromIndex = i;
        } else {

        }
      } else {
        peaks[peakNum] = chroms[i].l1Peaks[j];
        peaks[peakNum].chromIndex = i;

        peakNum++;
      }
    }
  }

  sort(peaks.begin(), peaks.end(), greaterFoldChange);

  for (i = 0; i < peakNum; i++) {
    if (peaks[i].isSignificant) {
      string chr(chroms[peaks[i].chromIndex].chromName);
      vector<uint32_t> tmp;
      tmp.push_back(peaks[i].peak);
      called_peak pk(peaks[i].start, peaks[i].end, peaks[i].qValue,
          peaks[i].qValue, peaks[i].l1Count, peaks[i].l2Count, tmp);
      insertPeak(chr, pk);
    }
  }

  return 1;
}

//PreProcessing: sort and get unique tags
void ccat::SortAndDedup(vector<chr_t>&chroms) {

  for (size_t i = 0; i < chroms.size(); i++) {
    ccat_aux::SortAndGetUniqueTags(chroms[i].l1PosTags);
    ccat_aux::SortAndGetUniqueTags(chroms[i].l2PosTags);
    ccat_aux::SortAndGetUniqueTags(chroms[i].l1NegTags);
    ccat_aux::SortAndGetUniqueTags(chroms[i].l2NegTags);
  }
}

void ccat::cmain(Reads& treads, Reads& creads, cmd_option_parser& option) {
  string projectName;
  vector<chr_t> chroms;
  size_t chromNum;
  double l1Ratio = 1;
  double l2Ratio = 1;
  size_t maxL1Count, maxL2Count;

  conf.fragmentSize=option.getExt_length();
  conf.bootstrapPass=option.bootstrapPass;
  conf.isStrandSensitiveMode = 0;
  conf.minCount=option.minCount;
  conf.minScore=option.minScore;
  conf.movingStep=option.movingStep;
  conf.outputNum=option.outputNum;
  conf.randomSeed=123456;
  conf.slidingWinSize=option.slidingWinSize;
  noise.setFragmentSize(conf.fragmentSize);
  config_validation(conf);

  projectName = option.getOutput_file();
  LoadData(chroms, treads, creads, chromNum);
  SortAndDedup(chroms);
  srand(conf.randomSeed);
  noise.ComputeNoiseRate(chroms, chromNum, l1Ratio, l2Ratio);
  pkFinder.PeakFinding(chroms, chromNum, l1Ratio, l2Ratio, maxL1Count,
      maxL2Count, conf);
  fdr.SignificanceAnalysis(chroms, chromNum, l1Ratio, l2Ratio, maxL1Count,
      maxL2Count, conf);
  UploadPeaks(chroms, chromNum, projectName.c_str(), conf);
}

void ccat::detectSummits(Reads& treatment_reads, Reads& control_reads,
    cmd_option_parser& option) {
  cmain(treatment_reads, control_reads, option);
}
