/*
 * ReadPairAux.h
 *
 *  Created on: May 6, 2012
 *      Author: xin
 */

#ifndef READPAIRAUX_H_
#define READPAIRAUX_H_
#include "common/stl_header.h"
#include "ReadPair.h"
#include "Read.h"
#include "short_reads/BlockedRead.h"
namespace reads {

class ReadPairPrinter {
public:
    void operator()(const ReadPair<Read>& rp);

};

bool R1IsInInChr(const ReadPair<Read>& rp, const std::string& chr);

bool R2IsInInChr(const ReadPair<Read>& rp, const std::string& chr);

bool eitherReadisInChr(const ReadPair<Read>& rp, const std::string& chr);

bool sameStartAndDir(ReadPair<BlockedRead>& rp1, ReadPair<BlockedRead>& rp2);

void filterSameStartAndDir(std::vector<ReadPair<BlockedRead> >&tags);

} /* namespace reads */
#endif /* READPAIRAUX_H_ */
