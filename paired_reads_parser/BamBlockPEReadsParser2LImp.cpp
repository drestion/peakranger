/*
 * BamBlockPEReadsParser2LImp.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/BamBlockPEReadsParser2LImp.h"
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
using namespace reads;

namespace parser {
namespace aux {

BamBlockPEReadsParser2LImp::BamBlockPEReadsParser2LImp() :
        mPRead(){
    mPRead.SetIsUnmapped(true);
    parserID = "BamBlockPEReadsParser2LImp";
}

BamBlockPEReadsParser2LImp::~BamBlockPEReadsParser2LImp() {

}

void BamBlockPEReadsParser2LImp::parse(const BamAlignment& read,
        PairEndedReads<BlockedRead>& rds, const RefVector& ref) {

    // when this is not added, unmapped reads will be in an empty
    // chr.
    //Will not enforce isGoodPERead()
    //that is not my job
    if (!(mPRead.IsMapped())) {
        mPRead = read;
    } else {
        if (mPRead.Position < read.Position) {
            insertRead(mPRead, read, rds, ref);
        } else {
            insertRead(read, mPRead, rds, ref);
        }
        mPRead.SetIsUnmapped(true);
    }
}

void BamBlockPEReadsParser2LImp::flush(
        reads::PairEndedReads<reads::BlockedRead>& reads,
        const BamTools::RefVector& ref) {
    //One reason you reach here:
    //the last read is one of a pair but it is unmapped
    //Another reason:
    //The file contains single end reads
    if (mPRead.IsMapped()) {
        //todo is this a good dummy read?
        BamAlignment bam;
        bam.Position = -1; // this is not good..
        bam.RefID = -1; // acceptable at most
        insertRead(mPRead, bam, reads, ref);
    }
}

} /* namespace aux */
} /* namespace parser */
