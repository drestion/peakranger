/*
 * StockBamMultipleDatasetsApp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef STOCKBAMMULTIPLEDATASETSAPP_H_
#define STOCKBAMMULTIPLEDATASETSAPP_H_

#include "short_reads/PairEndedReads.h"
#include "short_reads/BlockedRead.h"
#include "StockBamMultipleDatasetsAppImp.h"
#include "paired_reads_parser/BamBlockPEReadsParserImp.h"
#include "paired_reads_parser/BamBlockPEReadsParser.h"
#include "utils/Tracer.h"
namespace bam_app {

/*
 * Stock reads from multile bam files and
 * then analyze them as a whole. Cant do this
 * online, must be stock.
 */
class StockBamMultipleDatasetsApp {
public:
    StockBamMultipleDatasetsApp(utils::TimeStampTracer& tracer);
    virtual ~StockBamMultipleDatasetsApp();

    virtual void run(const std::vector<std::string>& files,std::ostream& os);

    virtual void parseReads(const std::vector<std::string>& files);
    virtual void report(std::ostream& os);
    void setImp(aux::StockBamMultipleDatasetsAppImp* imp);
    void setParserImp(parser::aux::BamBlockPEReadsParserImp* parser);
    std::string getAppId() const;
    void setAppId(std::string appId);

protected:
    utils::TimeStampTracer& tracer;
    std::string mAppID;

private:
    aux::StockBamMultipleDatasetsAppImp* mImp;
    parser::BamBlockPEReadsParser* mParser;
    std::vector<reads::PairEndedReads<reads::BlockedRead> > mReads;


};

} /* namespace bam_app */
#endif /* STOCKBAMMULTIPLEDATASETSAPP_H_ */
