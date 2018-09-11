/*
 * BarHist.cpp
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#include "BarHist.h"
#include "BarAux.h"
#include "short_reads/BlockedBamHit.h"

using namespace boost::icl;
using namespace ranger::concepts;
using namespace std;
using namespace reads;

namespace ranger {
namespace bar {

BarHist::BarHist() {
}

BarHist::~BarHist() {
}

void BarHist::getCount(const vector<Read>& hits, vector<RegionCount>& res) {
    accHits(hits, mCounter);
    ranger::bar::getCount(res, mCounter);
}

Read BarHist::extendReadToLength(const Read& r, int32_t length) {
    int32_t end = r.getEnd();
    if (length > 0) {
        end += (length - (end - r.getStart()));
    }
    end = end > 0 ? end : 0;
    return Read(r.getStart(), end, r.getChr().c_str(), r.getStrand());
}

void BarHist::add(const Read& read) {
    accHit(read, mCounter);
}

void BarHist::add(const std::vector<Read>& hits) {
    foreach(const Read& r, hits) {
        accHit(r, mCounter);
    }
}

void BarHist::add(const BamTools::BamAlignment& bam) {
    BamTools::RefVector foo;
    BlockedBamHit hit(bam, foo);
    vector<Read> rds;
    hit.getBlocks(rds);
    add(rds);
}

void operator +=(BarHist& b1, const BarHist& b2) {
    b1.mCounter += b2.mCounter;
}

void operator +=(std::map<std::string, BarHist>& b1,
        std::map<std::string, BarHist>& b2) {
    std::map<string, BarHist>::const_iterator it;
    string chr;
    for (it = b2.begin(); it != b2.end(); it++) {
        chr = it->first;
        b1[chr] += b2[chr];
    }
}

void BarHist::reset() {
    mCounter.clear();
}

} /* namespace bar */
} /* namespace ranger */
