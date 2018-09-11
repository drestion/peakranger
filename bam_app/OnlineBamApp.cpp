/*
 * OnlineBamApp.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "bam_app/OnlineBamApp.h"
#include "common/stl_header.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "utils/Guarded.h"
#include "common/ranger_debug.h"
using namespace BamTools;
using namespace std;
namespace bam_app {

OnlineBamApp::OnlineBamApp() :
        mImp(0), mCntToReport(10000000) {
}

OnlineBamApp::~OnlineBamApp() {
}

OnlineBamApp::OnlineBamApp(aux::OnlineBamAppImp* imp) :
        mImp(imp), mCntToReport(10000000) {
}

void OnlineBamApp::processReads(const std::string& file, ostream& os) {
    using namespace BamTools;
    using namespace utils;
    BamReader bam;
    BamAlignment read, mread;
    Guarded<FileNotGood> g(!(bam.Open(file)), file.c_str());
    const RefVector refvec = bam.GetReferenceData();
    while (bam.GetNextAlignment(read)) {
        mImp->process(read, refvec);
    }

    report(os);
}

void OnlineBamApp::report(std::ostream& os) {
    mImp->report(os);
}

} /* namespace bam_app */
