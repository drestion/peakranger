/*
 * StockBamMultipleDatasetsApp.cpp
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#include "bam_app/StockBamMultipleDatasetsApp.h"
using namespace reads;

namespace bam_app {

StockBamMultipleDatasetsApp::StockBamMultipleDatasetsApp(
        utils::TimeStampTracer& tra) :
        tracer(tra), mAppID("StockBamMultipleDatasetsApp"), mImp(0), mParser(0), mReads() {

}

StockBamMultipleDatasetsApp::~StockBamMultipleDatasetsApp() {
    if (mParser) {
        delete mParser;
    }
}

void StockBamMultipleDatasetsApp::run(const std::vector<std::string>& files,
        std::ostream& os) {
    parseReads(files);
    report(os);
}

void StockBamMultipleDatasetsApp::parseReads(
        const std::vector<std::string>& files) {
    foreach(std::string file, files) {
        tracer << "[" << mAppID << "]" << "Parsing " << file << "\n";
        mReads.push_back(PairEndedReads<BlockedRead>());
        mParser->parse(mReads.back(), file);
        tracer << "[" << mAppID << "]" << "Parsing Complete\n";
    }
}

void StockBamMultipleDatasetsApp::report(std::ostream& os) {
    tracer << "[" << mAppID << "]" << "Running app: " << mImp->getAppId();
    tracer << "\n";
    mImp->report(mReads, os);
    tracer << "[" << mAppID << "]" << "App complete\n";
}

void StockBamMultipleDatasetsApp::setImp(
        aux::StockBamMultipleDatasetsAppImp* imp) {
    mImp = imp;
}

void StockBamMultipleDatasetsApp::setParserImp(
        parser::aux::BamBlockPEReadsParserImp* parser) {
    mParser = new parser::BamBlockPEReadsParser(tracer, parser);
}

std::string StockBamMultipleDatasetsApp::getAppId() const {
    return mAppID;
}

void StockBamMultipleDatasetsApp::setAppId(std::string appId) {
    mAppID = appId;
}

} /* namespace bam_app */
