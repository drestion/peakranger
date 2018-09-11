/*
 * BamPEReadsParser.h
 *
 *  Created on: May 6, 2012
 *      Author: xin
 */

#ifndef BAMPEREADSPARSER_H_
#define BAMPEREADSPARSER_H_

#include "common/stl_header.h"
#include "short_reads/PairEndedReads.h"
#include "bamtools/BamAux.h"
#include "short_reads/Read.h"
namespace parser {

class BamPEReadsParser {
public:

    void parse(reads::PairEndedReads<reads::Read>& reads,
            const std::string& file);

private:

    void insertRead(const BamTools::BamAlignment& read,
            reads::PairEndedReads<reads::Read>& reads,
            const BamTools::RefVector& ref);

};

} /* namespace parser */
#endif /* BAMPEREADSPARSER_H_ */
