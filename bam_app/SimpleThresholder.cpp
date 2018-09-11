/*
 * SimpleThresholder.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "bam_app/SimpleThresholder.h"
#include "bar/BarHist.h"
#include "short_reads/PairEndedReadsAux.h"
#include "common/stl_header.h"
#include "utils/Tracer.h"
using namespace reads;
using namespace ranger::bar;
using namespace std;

namespace bam_app {
namespace aux {

SimpleThresholder::SimpleThresholder(utils::TimeStampTracer& tracer) :
        StockBamMultipleDatasetsAppImp(tracer), mLowerCutoff(0), mUpperCufoff(5) {
    setAppId("SimpleTresholder");
}

SimpleThresholder::~SimpleThresholder() {

}

inline bool SimpleThresholder::isInCutRegion(RegionCount & r, int32_t lhs,
        int32_t rhs) {
    return lhs < r.getCnt() && rhs > r.getCnt();
}

void SimpleThresholder::report(
        std::vector<PairEndedReads<BlockedRead> >& rdsVec, std::ostream& os) {

    if (rdsVec.size()) {
        PairEndedReads<BlockedRead>& rds = rdsVec[0];
        BarHist hist;
        vector<RegionCount> res;
        vector<Read> rrr;

        vector<string> chrs = rds.getChrs();

        foreach(string chr, chrs) {
            reads::getReadsInChr(rds, chr, rrr);
            tracer << "Got " << rrr.size() << " reads in " << chr << "\n";
            hist.getCount(rrr, res);
            foreach(RegionCount & r, res) {
                if (isInCutRegion(r, mLowerCutoff, mUpperCufoff)) {
                    os << chr << "\t" << r << "\n";
                }

            }
            rrr.resize(0);
            res.resize(0);
            hist.reset();
        }

    }
}

int32_t SimpleThresholder::getLowerCutoff() const {
    return mLowerCutoff;
}

void SimpleThresholder::setLowerCutoff(int32_t lowerCutoff) {
    mLowerCutoff = lowerCutoff;
}

int32_t SimpleThresholder::getUpperCufoff() const {
    return mUpperCufoff;
}

void SimpleThresholder::setUpperCufoff(int32_t upperCufoff) {
    mUpperCufoff = upperCufoff;
}

}
/* namespace aux */
} /* namespace bam_app */
