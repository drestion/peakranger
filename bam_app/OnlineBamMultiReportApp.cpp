/*
 * OnlineBamMultiReportApp.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: xfeng
 */

#include "bam_app/OnlineBamMultiReportApp.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "utils/Guarded.h"
#include "common/ranger_debug.h"
#include <iostream>
using namespace std;
namespace bam_app {

OnlineBamMultiReportApp::OnlineBamMultiReportApp(utils::TimeStampTracer& tracer) :
        tracer(tracer), mCnt(0), mCntToReport(10000000), mImp(0) {

}

OnlineBamMultiReportApp::~OnlineBamMultiReportApp() {

}

void OnlineBamMultiReportApp::report(std::string file_prefix) {
    tracer << "[" << "OnlineBamMultiReportApp" << "]" << "Running app: "
            << mImp->getAppId();
    tracer << "\n";
    mImp->report(file_prefix);
    tracer << "[" << "OnlineBamMultiReportApp" << "]" << "App complete\n";
}

OnlineBamMultiReportApp::OnlineBamMultiReportApp(utils::TimeStampTracer& tracer,
        aux::OnlineBamMultiReportAppImp* imp) :
        tracer(tracer), mCnt(0), mCntToReport(10000000), mImp(imp) {

}

void OnlineBamMultiReportApp::processReads(const std::string& file,
        const std::string& result_prefix) {
    using namespace BamTools;
    using namespace utils;
    BamReader bam;
    BamAlignment read, mread;
    Guarded<FileNotGood> g(!(bam.Open(file)), file.c_str());
    const RefVector refvec = bam.GetReferenceData();
    tracer << "[" << "OnlineBamMultiReportApp" << "]" << "App ID:"
            << mImp->getAppId() << "\n";
    while (bam.GetNextAlignment(read)) {
        mImp->process(read, refvec);
        if (++mCnt % mCntToReport == 0) {
            tracer << "[" << "OnlineBamMultiReportApp" << "]"
                    << "Reads processed:\t" << mCnt / 1000000 << "\tmillion\n";
        }
    }
    tracer << "[" << "OnlineBamMultiReportApp" << "]" << "Total reads:\t" << mCnt
            << "\n";
    report(result_prefix);
}

} /* namespace bam_app */
