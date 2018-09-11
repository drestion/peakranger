/*
 * StockBamMultipleDatasetsAppImp.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "bam_app/StockBamMultipleDatasetsAppImp.h"

namespace bam_app {
namespace aux {

StockBamMultipleDatasetsAppImp::StockBamMultipleDatasetsAppImp(utils::TimeStampTracer& tra) :
        tracer(tra),mAppID("StockBamMultiDatasetsAppImp") {

}

StockBamMultipleDatasetsAppImp::~StockBamMultipleDatasetsAppImp() {

}

void StockBamMultipleDatasetsAppImp::report(
        std::vector<reads::PairEndedReads<reads::BlockedRead> >& rdsVec,
        std::ostream& os) {
}

std::string StockBamMultipleDatasetsAppImp::getAppId() const {
    return mAppID;
}

void StockBamMultipleDatasetsAppImp::setAppId(std::string appId) {
    mAppID = appId;
}

} /* namespace aux */
} /* namespace bam_app */
