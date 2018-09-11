/*
 * BarAux.h
 *
 *  Created on: May 25, 2012
 *      Author: tania
 */

#ifndef BARAUX_H_
#define BARAUX_H_
#include "BarCounter.h"
#include "short_reads/Read.h"
#include "RegionCount.h"
namespace ranger {
namespace bar {

void accHits(const std::vector<reads::Read>& hits,
        boost::icl::BarCounter& counter);

void accHit(const reads::Read& hit, boost::icl::BarCounter& counter);
void getCount(std::vector<RegionCount>& res, boost::icl::BarCounter& counter);

void printBedGraphTrackline(std::ostream& os, bool isNewLine = true);

void printBedGraphTrackline(std::ostream& of, const char* _name = "",
        std::vector<uint32_t> col = std::vector<uint32_t>(0), bool isNewLine =
                true);

}
}
#endif /* BARAUX_H_ */
