/*
 * PairEndedReadsAux.h
 *
 *  Created on: May 27, 2012
 *      Author: tania
 */

#ifndef PAIRENDEDREADSAUX_H_
#define PAIRENDEDREADSAUX_H_

#include "PairEndedReads.h"
#include "BlockedRead.h"
#include "Read.h"
#include "common/stl_header.h"
namespace reads {

void getReadsInChr(PairEndedReads<BlockedRead>& rds, const std::string& chr,
		std::vector<Read>& results);
}
#endif /* PAIRENDEDREADSAUX_H_ */
