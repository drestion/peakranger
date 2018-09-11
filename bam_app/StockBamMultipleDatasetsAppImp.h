/*
 * StockBamMultipleDatasetsAppImp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef STOCKBAMMULTIPLEDATASETSAPPIMP_H_
#define STOCKBAMMULTIPLEDATASETSAPPIMP_H_
#include "short_reads/PairEndedReads.h"
#include "short_reads/BlockedRead.h"
#include "utils/Tracer.h"
#include <string>
namespace bam_app {
namespace aux {

class StockBamMultipleDatasetsAppImp {
public:
    StockBamMultipleDatasetsAppImp(utils::TimeStampTracer& tracer);
    virtual ~StockBamMultipleDatasetsAppImp();

    virtual
    void report(std::vector<reads::PairEndedReads<reads::BlockedRead> >& rdsVec,
            std::ostream& os);
    std::string getAppId() const;
    void setAppId(std::string appId);

protected:
    utils::TimeStampTracer& tracer;
    std::string mAppID;
};

} /* namespace aux */
} /* namespace bam_app */
#endif /* STOCKBAMMULTIPLEDATASETSAPPIMP_H_ */
