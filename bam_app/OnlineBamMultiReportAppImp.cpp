/*
 * OnlineBamMultiReportAppImp.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: xfeng
 */

#include "bam_app/OnlineBamMultiReportAppImp.h"

namespace bam_app {
namespace aux {
OnlineBamMultiReportAppImp::OnlineBamMultiReportAppImp(
        utils::TimeStampTracer& tracer) :
        mAppID("OnlineBamMultiReportAppImp"), tracer(tracer) {

}

OnlineBamMultiReportAppImp::~OnlineBamMultiReportAppImp() {

}

void OnlineBamMultiReportAppImp::process(const BamTools::BamAlignment& read,
        const BamTools::RefVector& ref) {
}

void OnlineBamMultiReportAppImp::report(std::string& file_prefix) {
}

std::string OnlineBamMultiReportAppImp::getAppId() const {
    return mAppID;
}

void OnlineBamMultiReportAppImp::setAppId(std::string appId) {
    mAppID = appId;
}
}
} /* namespace bam_app */
