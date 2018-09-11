/*
 * BamBlockPEReadsParser.h
 *
 *  Created on: May 10, 2012
 *      Author: xin
 */

#ifndef BAMBLOCKPEREADSPARSER_H_
#define BAMBLOCKPEREADSPARSER_H_

#include "short_reads/PairEndedReads.h"
#include "short_reads/BlockedRead.h"
#include "common/stl_header.h"
#include "BamBlockPEReadsParserImp.h"
#include "utils/Tracer.h"
namespace parser {

class BamBlockPEReadsParser {
public:
    BamBlockPEReadsParser(utils::TimeStampTracer& tracer);
    ~BamBlockPEReadsParser();
    BamBlockPEReadsParser(utils::TimeStampTracer& tracer,aux::BamBlockPEReadsParserImp* parser);

    void parse(reads::PairEndedReads<reads::BlockedRead>& reads,
            const std::string& file);

    aux::BamBlockPEReadsParserImp* getParser() const {
        return mParser;
    }

    void setParser(aux::BamBlockPEReadsParserImp* parser) {
        mParser = parser;
    }

protected:
    utils::TimeStampTracer& tracer;
private:
    aux::BamBlockPEReadsParserImp* mParser;
    uint32_t mCnt;
};

} /* namespace parser */
#endif /* BAMBLOCKPEREADSPARSER_H_ */
