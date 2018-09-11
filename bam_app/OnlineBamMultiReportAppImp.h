/*
 * OnlineBamMultiReportAppImp.h
 *
 *  Created on: Jun 20, 2012
 *      Author: xfeng
 */

#ifndef ONLINEBAMMULTIREPORTAPPIMP_H_
#define ONLINEBAMMULTIREPORTAPPIMP_H_
#include <string>
#include <stdint.h>
#include "bamtools/BamAux.h"
#include "utils/Tracer.h"
namespace bam_app {
namespace aux {
class OnlineBamMultiReportAppImp {
public:
    OnlineBamMultiReportAppImp(utils::TimeStampTracer& tracer);
    virtual ~OnlineBamMultiReportAppImp();

    virtual void process(const BamTools::BamAlignment & read,
            const BamTools::RefVector & ref);
    virtual void report(std::string& file_prefix);
    std::string getAppId() const;
    void setAppId(std::string appId);

protected:
    std::string mAppID;
    utils::TimeStampTracer& tracer;
};
} /* namesapce aux */
} /* namespace bam_app */
#endif /* ONLINEBAMMULTIREPORTAPPIMP_H_ */
