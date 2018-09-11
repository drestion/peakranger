/*
 * CCATNoiseRate.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#include "CCATNoiseRate.h"
#include "common/ranger_debug.h"
using namespace std;
namespace ccat_aux {
const size_t CCATNoiseRate::BINSIZE = 1000;
const size_t CCATNoiseRate::MAXITERATION = 20;
const size_t CCATNoiseRate::MINCOUNT = 10000;

CCATNoiseRate::CCATNoiseRate() :
        fragSize(0) {

}

CCATNoiseRate::~CCATNoiseRate() {

}

void CCATNoiseRate::setFragmentSize(const size_t& fs) {
    fragSize = fs;
}

double CCATNoiseRate::ComputeNoiseRate(
        const std::vector<ccat_aux::chr_t> & chroms, size_t chromNum,
        double & l1Ratio, double & l2Ratio) const {
    LOG_DEBUG1("Entering CCATNoiseRate::ComputeNoiseRate");
    size_t i;
    double nr, prevNr;
    double sum;
    int count, flag;
    int l1TagNum = 0, l2TagNum = 0;

    nr = 1.0;
    prevNr = 1.0;

    sum = 0;
    count = 0;
    flag = 0;

    for (i = 0; i < MAXITERATION; i++) {
        LOG_DEBUG2("In iteration " << i);
        ComputeNoiseRateInOneIteration(chroms, chromNum, nr);
        if (!(nr > 0 && nr < 1)) {
            nr = 1;
            break;
        }
        if (nr > prevNr) {
            flag = 1;
        }

        prevNr = nr;

        if (flag) {
            sum += nr;
            count++;
        }

    }

    if (count > 0) {
        nr = sum / count;
    }

    for (i = 0; i < chromNum; i++) {
        l1TagNum += chroms[i].l1PosTags.size() + chroms[i].l1NegTags.size();
        l2TagNum += chroms[i].l2PosTags.size() + chroms[i].l2NegTags.size();
    }

    ComputeSampleRatio(nr, l1TagNum, l2TagNum, l1Ratio, l2Ratio);
    LOG_DEBUG1("Leaving CCATNoiseRate::ComputeNoiseRate");
    return nr;
}

void CCATNoiseRate::ComputeNoiseRateInOneIteration(
        const std::vector<ccat_aux::chr_t> & chroms, size_t chromNum,
        double & noiseRate) const {
    MARK_FUN("CCATNoiseRate::ComputeNoiseRateInOneIteration");

    LOG_DEBUG2("chromNum "<<chromNum);

    LOG_DEBUG2("noiseRate "<<noiseRate);LOG_DEBUG2("MINCOUNT "<<MINCOUNT);
    size_t i, j;
    size_t l1TagNum = 0, l2TagNum = 0;
    double l1Ratio, l2Ratio;
    int thresh, tmpL1Count, tmpL2Count, sum1, sum2;
    vector<vector<bin_t> > bins(chromNum);
    vector<size_t> binNum(chromNum, 0);

    for (i = 0; i < chromNum; i++) {
        l1TagNum += chroms[i].l1PosTags.size() + chroms[i].l1NegTags.size();
        l2TagNum += chroms[i].l2PosTags.size() + chroms[i].l2NegTags.size();
    }

    LOG_DEBUG2("l1TagNum "<<l1TagNum);

    LOG_DEBUG2("l2TagNum "<<l2TagNum);

    //+1 to make sure indexing is legal.
    vector<size_t> hist(l2TagNum + 1, 0);
    ComputeSampleRatio(noiseRate, l1TagNum, l2TagNum, l1Ratio, l2Ratio);

    LOG_DEBUG2("l1Ratio "<<l1Ratio);LOG_DEBUG2("l2Ratio "<<l2Ratio);

    AllocBinMem(chroms, chromNum, bins, binNum);

    for (i = 0; i < chromNum; i++) {
        LOG_DEBUG3("In chr number "<<i);
        CountFragmentsInOneChrom(chroms.at(i), bins[i], binNum[i], l1Ratio,
                l2Ratio, fragSize);
        for (j = 0; j < binNum[i]; j++) {
            tmpL1Count = bins[i][j].l1Counts;
            tmpL2Count = bins[i][j].l2Counts;
            if (tmpL1Count || tmpL2Count) {
                LOG_DEBUG4("Using bin counts from "<<j);

                LOG_DEBUG4("tmpL1Count: "<<tmpL1Count);LOG_DEBUG4("tmpL2Count: "<<tmpL2Count);

                LOG_DEBUG4("Updated hist[ "<<(tmpL1Count < tmpL2Count ? tmpL2Count : tmpL1Count)<<"]");

                LOG_DEBUG4("        from "<<hist[(tmpL1Count < tmpL2Count ? tmpL2Count : tmpL1Count)]);
                hist[tmpL1Count < tmpL2Count ? tmpL2Count : tmpL1Count] +=
                        tmpL2Count;
                LOG_DEBUG4("        to   "<<hist[(tmpL1Count < tmpL2Count ? tmpL2Count : tmpL1Count)]);
            }
        }
    }

    sum1 = 0;

    for (i = 0; i < l2TagNum; i++) {
        sum1 += hist[i];

        if (((double) sum1 > MINCOUNT)
                && ((double) sum1 > l2TagNum * l2Ratio * 0.5)) {
            break;
        }

    }

    LOG_DEBUG3(" before sum1 : "<<sum1);

    if (i == l2TagNum) {

        return;
    }

    thresh = i;

    LOG_DEBUG3("thresh "<<thresh);

    sum1 = 0;
    sum2 = 0;

    for (i = 0; i < chromNum; i++) {
        for (j = 0; j < binNum[i]; j++) {
            tmpL1Count = bins[i][j].l1Counts;
            tmpL2Count = bins[i][j].l2Counts;

            if ((tmpL1Count < tmpL2Count ? tmpL2Count : tmpL1Count) <= thresh) {
                sum1 += tmpL1Count;
                sum2 += tmpL2Count;
            }
        }
    }

    noiseRate *= (double) sum1 / sum2;
    LOG_DEBUG3("sum1 : "<<sum1);LOG_DEBUG3("sum2 : "<<sum2);

    LOG_DEBUG3("Final noiseRate "<<noiseRate);

}

void CCATNoiseRate::CountFragmentsInOneChrom(const ccat_aux::chr_t & chrom,
        std::vector<ccat_aux::bin_t> & bins, int binNum, double l1Ratio,
        double l2Ratio, size_t fragmentSize) const {
    MARK_FUN("CCATNoiseRate::CountFragmentsInOneChrom");

    LOG_DEBUG2("fragmentSize: "<<fragmentSize);LOG_DEBUG2("l1Ratio: "<<l1Ratio);

    LOG_DEBUG2("l2Ratio: "<<l2Ratio);LOG_DEBUG2("BINSIZE: "<<BINSIZE);

    size_t j;
    int tmpIndex;

    //count fragments for L1
    for (j = 0; j < binNum; j++) {
        bins[j].l1Counts = 0;
        bins[j].l2Counts = 0;
    }

    for (j = 0; j < chrom.l1PosTags.size(); j++) {
        if (rand() > RAND_MAX * l1Ratio) {
            continue;
        }

        tmpIndex = (chrom.l1PosTags[j] + fragmentSize / 2) / BINSIZE;

        tmpIndex = tmpIndex >= binNum ? binNum - 1 : tmpIndex;
        tmpIndex = tmpIndex < 0 ? 0 : tmpIndex;

        LOG_DEBUG5("Got pos tag: "<<chrom.l1PosTags[j]<<" and ind: "<<tmpIndex);
        bins[tmpIndex].l1Counts++;
    }

    for (j = 0; j < chrom.l1NegTags.size(); j++) {
        if (rand() > RAND_MAX * l1Ratio) {
            continue;
        }

        tmpIndex = (chrom.l1NegTags[j] - fragmentSize / 2) / BINSIZE;

        tmpIndex = tmpIndex >= binNum ? binNum - 1 : tmpIndex;
        tmpIndex = tmpIndex < 0 ? 0 : tmpIndex;
        LOG_DEBUG5("Got neg tag: "<<chrom.l1NegTags[j]<<" and ind: "<<tmpIndex);
        bins[tmpIndex].l1Counts++;
    }

    for (j = 0; j < chrom.l2PosTags.size(); j++) {
        if (rand() > RAND_MAX * l2Ratio) {
            continue;
        }

        tmpIndex = (chrom.l2PosTags[j] + fragmentSize / 2) / BINSIZE;

        tmpIndex = tmpIndex >= binNum ? binNum - 1 : tmpIndex;
        tmpIndex = tmpIndex < 0 ? 0 : tmpIndex;

        bins[tmpIndex].l2Counts++;
    }

    for (j = 0; j < chrom.l2NegTags.size(); j++) {
        if (rand() > RAND_MAX * l2Ratio) {
            continue;
        }

        tmpIndex = (chrom.l2NegTags[j] - fragmentSize / 2) / BINSIZE;

        tmpIndex = tmpIndex >= binNum ? binNum - 1 : tmpIndex;
        tmpIndex = tmpIndex < 0 ? 0 : tmpIndex;

        bins[tmpIndex].l2Counts++;
    }

}
void CCATNoiseRate::ComputeSampleRatio(double noiseRate, int l1TagNum,
        int l2TagNum, double& l1Ratio, double& l2Ratio) const {
    MARK_FUN("CCATNoiseRate::ComputeSampleRatio");LOG_DEBUG4(" noiseRate:"<<noiseRate);
    if (l1TagNum * noiseRate > l2TagNum) {
        l1Ratio = (double) (l2TagNum) / (l1TagNum * noiseRate);
        l2Ratio = 1.0;
    } else {
        l2Ratio = (double) (l1TagNum * noiseRate) / l2TagNum;
        l1Ratio = 1.0;

    }
}
void CCATNoiseRate::AllocBinMem(const vector<chr_t>&chroms, int chromNum,
        std::vector<std::vector<bin_t> >& bins,
        std::vector<size_t>& binNum) const {
    MARK_FUN("Entering CCATNoiseRate::AllocBinMem");
    int i;

    for (i = 0; i < chromNum; i++) {
        if (chroms[i].chromSize > 0) {
            LOG_DEBUG3("IN AllocBinMem: at chr "<<i<<" chromSize : "<<chroms[i].chromSize);
            binNum[i] = chroms[i].chromSize / BINSIZE + 1;
            LOG_DEBUG3("IN AllocBinMem: at chr "<<i<<" num_of_bin : "<<chroms[i].chromSize / BINSIZE + 1);
            vector<bin_t> tmp(binNum[i]);
            bins[i].swap(tmp);
        }
    }

}
} /* namespace ccat_aux */
