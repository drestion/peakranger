/*
 * BamBlockPEReadsParserImp.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/BamBlockPEReadsParserImp.h"

namespace parser {
namespace aux {

BamBlockPEReadsParserImp::BamBlockPEReadsParserImp():parserID("Null Imp") {
}

BamBlockPEReadsParserImp::~BamBlockPEReadsParserImp() {
}

void BamBlockPEReadsParserImp::parse(const BamTools::BamAlignment& read,
        reads::PairEndedReads<reads::BlockedRead>& reads,
        const BamTools::RefVector& ref) {
}

void BamBlockPEReadsParserImp::flush(
        reads::PairEndedReads<reads::BlockedRead>& reads,
        const BamTools::RefVector& ref) {
}

} /* namespace aux */
} /* namespace parser */
