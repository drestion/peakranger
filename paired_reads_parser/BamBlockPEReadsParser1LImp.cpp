/*
 * BamBlockPEReadsParser1LImp.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/BamBlockPEReadsParser1LImp.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "short_reads/BlockedBamFormatAux.h"
#include "short_reads/bamFormatAux.h"
#include "short_reads/BamHit.h"
#include "short_reads/ReadPair.h"
#include "short_reads/BlockedBamHit.h"
#include "BamParserAux.h"
#include "common/ranger_debug.h"
using namespace BamTools;
using namespace std;
using namespace reads;
namespace parser {
namespace aux {

BamBlockPEReadsParser1LImp::BamBlockPEReadsParser1LImp()
       {
    mPRead.RefID = -1;
    parserID="BamBlockPEReadsParser1LImp";
}

BamBlockPEReadsParser1LImp::~BamBlockPEReadsParser1LImp() {
}

void BamBlockPEReadsParser1LImp::parse(const BamTools::BamAlignment& read,
        reads::PairEndedReads<reads::BlockedRead>& rds,
        const BamTools::RefVector& ref) {
//	cout <<"[BamBlockPEReadsParser1LImp]parse "<<read.Position<<"\n";
    insertRead(read, mPRead, rds, ref);
}

void BamBlockPEReadsParser1LImp::flush(
        reads::PairEndedReads<reads::BlockedRead>& reads,
        const BamTools::RefVector& ref) {
}


} /* namespace aux */
} /* namespace parser */
