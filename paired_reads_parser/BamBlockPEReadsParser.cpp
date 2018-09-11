/*
 * BamBlockPEReadsParser.cpp
 *
 *  Created on: May 10, 2012
 *      Author: xin
 */

#include "BamBlockPEReadsParser.h"
#include "common/ranger_debug.h"
#include "common/boost_header.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "utils/Guarded.h"
#include "BamBlockPEReadsParserhotelImp.h"
using namespace reads;
using namespace std;
namespace parser {

BamBlockPEReadsParser::BamBlockPEReadsParser(utils::TimeStampTracer& tracer) :
        tracer(tracer), mParser(), mCnt(0) {
}

BamBlockPEReadsParser::BamBlockPEReadsParser(utils::TimeStampTracer& tracer,
        aux::BamBlockPEReadsParserImp* parser) :
        tracer(tracer), mParser(parser), mCnt(0) {
}

BamBlockPEReadsParser::~BamBlockPEReadsParser() {
}

void parser::BamBlockPEReadsParser::parse(PairEndedReads<BlockedRead>& reads,
        const std::string& file) {
    using namespace BamTools;
    using namespace utils;
    BamReader bam;
    BamAlignment read, mread;
    Guarded<FileNotGood> g(!(bam.Open(file)), file.c_str());
    const RefVector refvec = bam.GetReferenceData();
    tracer << "[" << "BamBlockPEReadsParser" << "]" << "Paser implementation:"
            << mParser->getParserId() << "\n";
    while (bam.GetNextAlignment(read)) {
        //todo: should be a general vector of filters of reads
        if (read.IsMapped()) {
            mParser->parse(read, reads, refvec);
        }
        if (++mCnt % 10000000 == 0) {
            tracer << "[" << "BamBlockPEReadsParser" << "]";
            tracer << "Reads processed:\t" << mCnt / 1000000 << "\tmillion\n";
        }
    }
    tracer << "[" << "BamBlockPEReadsParser" << "]";
    tracer << "Total reads:\t" << mCnt << "\n";
    mParser->flush(reads, refvec);
}
}

