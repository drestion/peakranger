/*
 * WigPEOnline.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "WigPEOnline.h"
#include <string>
#include "short_reads/bamFormatAux.h"
#include "short_reads/BlockedBamFormatAux.h"
#include "short_reads/Read.h"
#include "short_reads/BlockedBamHit.h"
#include "common/boost_header.h"
#include "bar/BarCounter.h"
#include "bar/BGHistStream.h"
#include "bar/BarAux.h"

using namespace boost::icl;
using namespace boost;
using namespace std;
using namespace ranger::bar;
using namespace reads;

namespace bam_app {
namespace aux {

WigPEOnline::WigPEOnline(utils::TimeStampTracer& tracer) :
        OnlineBamMultiReportAppImp(tracer), mCounter(), mNCounter(), mExt(0), mbSplitByStrand(
                false), mbSplitByChr(false), mPosCnt(0), mNegCnt(0), gz(false) {
    setAppId("WigPEOnline");
}

WigPEOnline::~WigPEOnline() {
}

void WigPEOnline::process(const BamTools::BamAlignment& read,
        const BamTools::RefVector& ref) {
    if (read.IsMapped()) {
        string chr = reads::getR1Chr(read, ref);
        vector<Read> rds;

        if (mExt) {
            extendBamRead(read, rds, mExt);
        } else {
            parseBamBlocks(rds, read, ref);
        }
        if (read.IsReverseStrand()) {
            mNegCnt++;
            mNCounter[chr].add(rds);
        } else {
            mPosCnt++;
            mCounter[chr].add(rds);
        }
    }
}

void WigPEOnline::report(string& _pre) {
    std::map<string, BarHist>::iterator it;
    if (isSplitByChr() && isSplitByStrand()) {
        tracer << "[" << getAppId() << "]" << "Split by Chr and Strand\n";
        tracer << "[" << getAppId() << "]" << "Exporting Pos reads...\n";
        outputAllByChrWithHeader(mCounter, _pre, string("_Pos.wig"), gz);
        tracer << "[" << getAppId() << "]" << "Complete\n";
        tracer << "[" << getAppId() << "]" << "Exporting Neg reads...\n";
        outputAllByChrWithHeader(mNCounter, _pre, string("_Neg.wig"), gz);
        tracer << "[" << getAppId() << "]" << "Complete\n";
    } else if (!isSplitByChr() && isSplitByStrand()) {
        tracer << "[" << getAppId() << "]" << "Split by Strand\n";
        tracer << "[" << getAppId() << "]" << "Exporting Pos reads...\n";
        outputAllWithHeader(mCounter, _pre + string("_Pos.wig"), gz);
        tracer << "[" << getAppId() << "]" << "Complete\n";
        tracer << "[" << getAppId() << "]" << "Exporting Neg reads...\n";
        outputAllWithHeader(mNCounter, _pre + string("_Neg.wig"), gz);
        tracer << "[" << getAppId() << "]" << "Complete\n";
    } else if (!isSplitByStrand()) {
        std::map<string, BarHist> all(mCounter);
        all += mNCounter;
        if (isSplitByChr()) {
            tracer << "[" << getAppId() << "]" << "Split by Chr\n";
            tracer << "[" << getAppId() << "]" << "Exporting all reads...\n";
            outputAllByChrWithHeader(all, _pre, string(".wig"), gz);
            tracer << "[" << getAppId() << "]" << "Complete\n";
        } else {
            tracer << "[" << getAppId() << "]" << "Exporting all reads...\n";
            outputAllWithHeader(all, _pre + string(".wig"), gz);
            tracer << "[" << getAppId() << "]" << "Complete\n";
        }
    }
}

} /* namespace aux */
} /* namespace bam_app */
