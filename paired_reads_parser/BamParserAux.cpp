/*
 * BamParserAux.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "BamParserAux.h"
#include "short_reads/BlockedBamFormatAux.h"
#include "utils/Guarded.h"
#include "bamtools/BamReader.h"
#include "common/ranger_debug.h"
using namespace BamTools;
using namespace utils;
using namespace std;
namespace parser {
namespace aux {

void insertRead(const BamAlignment& read, const BamAlignment& mread,
        reads::PairEndedReads<reads::BlockedRead>& rds, const RefVector& ref) {
    using namespace reads;
    ReadPair<BlockedRead> rp;
    rp.add(getBlockedRead1FromPEBamRead(read, ref),
            getBlockedRead1FromPEBamRead(mread, ref));
//	cout <<"[insertRead]after"<< getBlockedRead1FromPEBamRead(read, ref)<<"\n";
//    cout <<"and added:"<< getBlockedRead1FromPEBamRead(mread, ref)<<"\n";
//    cout <<"R1 chr "<<rp.r1().getChr()<<"\n";
//    cout <<"R2 chr "<<rp.r2().getChr()<<"\n";
    rds.addReadPair(rp);
//    cout <<"[insertRead]after addReadPair"<<"\n";
}

}
}

void parser::fetchLines(vector<BamAlignment>& result, uint32_t n,
        const std::string& file) {
    BamReader bam;
    BamAlignment read;
    Guarded<FileNotGood> g(!(bam.Open(file)), file.c_str());
    const RefVector refvec = bam.GetReferenceData();
    while (bam.GetNextAlignment(read) && n) {
        result.push_back(read);
//        cout << "read " << n << "\t" << read << "\n";
        n--;
    }
}

