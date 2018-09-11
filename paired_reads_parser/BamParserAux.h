/*
 * BamParserAux.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef BAMPARSERAUX_H_
#define BAMPARSERAUX_H_
#include "bamtools/BamAux.h"
#include "short_reads/BlockedRead.h"
#include "short_reads/PairEndedReads.h"
#include "short_reads/ReadPair.h"
#include "common/stl_header.h"
namespace parser {
namespace aux {

void insertRead(const BamTools::BamAlignment& read,
        const BamTools::BamAlignment& mread,
        reads::PairEndedReads<reads::BlockedRead>& rds,
        const BamTools::RefVector& ref);
}

void fetchLines(std::vector<BamTools::BamAlignment>& result, uint32_t n,const std::string& file);
}

#endif /* BAMPARSERAUX_H_ */
