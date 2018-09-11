/*
 * BlockedBamHit.cpp
 *
 *  Created on: May 10, 2012
 *      Author: xin
 */

#include "BlockedBamHit.h"
#include "CigarString.h"
#include "Read.h"
#include "Strand.h"
#include "BlockedBamFormatAux.h"
#include "bamFormatAux.h"
#include "concepts/Region.h"
#include "CigarStringFactory.h"
#include "common/ranger_debug.h"
using namespace std;
using namespace BamTools;

namespace reads {

BlockedBamHit::BlockedBamHit() {
}

BlockedBamHit::~BlockedBamHit() {

}

BlockedBamHit::BlockedBamHit(const BamAlignment & bam, const RefVector & ref) {
    parseBamBlocks(mBlocks, bam,ref);
}

} /* namespace reads */
