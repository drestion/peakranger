/*
 * OnlineBamMultiReportApp.h
 *
 *  Created on: Jun 20, 2012
 *      Author: xfeng
 */

#ifndef ONLINEBAMMULTIREPORTAPP_H_
#define ONLINEBAMMULTIREPORTAPP_H_
#include "bam_app/OnlineBamMultiReportAppImp.h"
#include "utils/Tracer.h"
namespace bam_app {

class OnlineBamMultiReportApp {
public:
    OnlineBamMultiReportApp(utils::TimeStampTracer& tracer);
    virtual ~OnlineBamMultiReportApp();
    OnlineBamMultiReportApp(utils::TimeStampTracer& tracer,aux::OnlineBamMultiReportAppImp* imp);
    virtual void report(std::string file_prefix);
    virtual void processReads(const std::string& data,
            const std::string& file_prefix);
protected:
    utils::TimeStampTracer& tracer;

private:
    aux::OnlineBamMultiReportAppImp* mImp;
    uint32_t mCnt;
    uint32_t mCntToReport;
};

} /* namespace bam_app */
#endif /* ONLINEBAMMULTIREPORTAPP_H_ */
