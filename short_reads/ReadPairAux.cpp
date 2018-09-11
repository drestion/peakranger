/*
 * ReadPairAux.cpp
 *
 *  Created on: May 6, 2012
 *      Author: xin
 */

#include "ReadPairAux.h"
#include <iostream>
using namespace std;
void reads::ReadPairPrinter::operator ()(const ReadPair<Read>& rp) {
    cout << rp.r1();
    cout << rp.r2();
}

namespace reads {

bool R1IsInInChr(const ReadPair<Read>& rp, const std::string& chr) {
    return rp.r1().getChr() == chr;
}

bool R2IsInInChr(const ReadPair<Read>& rp, const std::string& chr) {
    return rp.r2().getChr() == chr;
}

bool eitherReadisInChr(const ReadPair<Read>& rp, const std::string& chr) {
    return R1IsInInChr(rp, chr) || R2IsInInChr(rp, chr);
}

bool sameStartAndDir(ReadPair<BlockedRead>& rp1, ReadPair<BlockedRead>& rp2) {

    return rp1.r1().getDir() == rp2.r1().getDir()
            && rp1.r2().getDir() == rp2.r2().getDir()
            && rp1.r1().getStart() == rp2.r1().getStart()
            && rp1.r2().getStart() == rp2.r2().getStart();

}

void filterSameStartAndDir(std::vector<ReadPair<BlockedRead> > tags) {
    tags.resize(
            std::unique(tags.begin(), tags.end(), sameStartAndDir)
                    - tags.begin());
}

} /* namespace reads */
