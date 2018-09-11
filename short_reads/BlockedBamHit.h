/*
 * BlockedBamHit.h
 *
 *  Created on: May 10, 2012
 *      Author: xin
 */

#ifndef BLOCKEDBAMHIT_H_
#define BLOCKEDBAMHIT_H_

#include "BlockedRead.h"
#include "bamtools/BamAux.h"
namespace reads {

class BlockedBamHit: public BlockedRead {
public:
    BlockedBamHit();
    virtual ~BlockedBamHit();

    BlockedBamHit(const BamTools::BamAlignment& bam, const BamTools::RefVector& ref);

};

} /* namespace reads */
#endif /* BLOCKEDBAMHIT_H_ */
