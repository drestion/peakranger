/*
 * LibraryComplexity.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: xfeng
 */

#include "bam_app/LibraryComplexity.h"
#include "short_reads/ReadPairAux.h"
using namespace std;
using namespace reads;
namespace bam_app {
namespace aux {

LibraryComplexity::LibraryComplexity(utils::TimeStampTracer& tra) :
        StockBamMultipleDatasetsAppImp(tra), mIsPE(true) {
    setAppId("Library Complexity");
}

LibraryComplexity::~LibraryComplexity() {

}

void LibraryComplexity::report(vector<PairEndedReads<BlockedRead> >& rdsVec,
        ostream& os) {
    tracer << "[" << mAppID << "]" << "Start " << "\n";
    if (rdsVec.size()) {
        PairEndedReads<BlockedRead>& rds(rdsVec[0]);

        size_t be = 0;
        size_t af = 0;
        foreach(string chr, rds.getChrs()) {
            tracer << "[" << mAppID << "]" << "Processing " << chr << "\n";
            be += rds.endOf(chr) - rds.beginOf(chr);
            tracer << "[" << mAppID << "]" << "Total reads in " << chr << ":\t";
            tracer << rds.endOf(chr) - rds.beginOf(chr) << "\n";
            vector<ReadPair<BlockedRead> > ga;
            rds.getReadPairs(chr.c_str(), ga);
            ga.resize(
                    std::unique(ga.begin(), ga.end(), sameStartAndDir)
                            - ga.begin());
            tracer << "[" << mAppID << "]" << "Unique reads in " << chr << ":\t";
            tracer << ga.size() << "\n";
            af += ga.size();
        }

        if (rdsVec.size()) {
            if (mIsPE) {
                os << "Unique read pairs:        " << af << "\n";
                os << "Total read pairs:         " << be << "\n";
            } else {
                os << "Unique reads:             " << af << "\n";
                os << "Total reads:              " << be << "\n";
            }
            os << "Library complexity:       " << setprecision(3)
                    << af * 100.0 / be << "%\n";
        } else {
            os << "No reads after filtering.\n";
        }
    }
}

} /* namespace aux */
} /* namespace bam_app */
