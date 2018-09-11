/*
 * BamPEReadsParser.cpp
 *
 *  Created on: May 6, 2012
 *      Author: xin
 */

#include "BamPEReadsParser.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "short_reads/bamFormatAux.h"
#include "short_reads/BamHit.h"
#include "short_reads/ReadPair.h"

using namespace std;
using namespace BamTools;
using namespace reads;

namespace parser {

void insertPairedRead(const BamAlignment& read, PairEndedReads<Read>& reads, const BamTools::RefVector& ref) {
    BamHit bam1;
    BamHit bam2;
    getBamHitsFromPEBamRead(read, ref, bam1, bam2);
    ReadPair<Read> rp;
    rp.add(bam1, bam2);
    reads.addReadPair(rp);
}

void BamPEReadsParser::insertRead(const BamAlignment& read, PairEndedReads<Read>& reads, const BamTools::RefVector& ref) {
    if (isFirstPEGoodRead(read)) {
        insertPairedRead(read, reads, ref);
    } else if (isGoodSERead(read)) {
        throw NotPairEndBamRead("The bam file contains single end reads.");
    } else {
        processAbnormalRead(read, ref);
    }
}

void BamPEReadsParser::parse(PairEndedReads<Read>& reads, const std::string& filename) {
    //todo: This parser tries to get info for paired reads from the first mate
    //todo: , which is good except for getting the EXACT end location of the second
    //todo:  mate. To solve this, a better parser should be able pair reads with the
    //todo: the same read name and treat each of them as the first mate.
    BamReader bam;
    BamAlignment read;
    if (!(bam.Open(filename))) {
        throw FileNotGood(filename.c_str());
    }
    const BamTools::RefVector refvec = bam.GetReferenceData();
    while (bam.GetNextAlignment(read)) {
        insertRead(read, reads, refvec);
    }
}

}
/* namespace parser */
