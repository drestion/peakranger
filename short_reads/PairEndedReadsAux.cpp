/*
 * PairEndedReadsAux.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */


#include "PairEndedReadsAux.h"

namespace reads{
void getReadsInChr(PairEndedReads<BlockedRead>& rds, const std::string& chr,
        std::vector<Read>& results) {
    using namespace std;

    std::vector<ReadPair<BlockedRead> >::iterator l, r;
    vector<string> chrs = rds.getChrs();
    foreach( string& c, chrs) {

        l = rds.beginOf(c);
        r = rds.endOf(c);
        while (l != r) {
            ReadPair<BlockedRead> rp = *l++;
            if (rp.r1().getChr() == chr) {
                rp.r1().getBlocks(results);
            }
            if (rp.r2().getChr() == chr) {
                rp.r2().getBlocks(results);
            }
        }

    }

}

}

