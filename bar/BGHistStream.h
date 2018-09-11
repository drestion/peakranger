/*
 * BGHistStream.h
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#ifndef BGHISTSTREAM_H_
#define BGHISTSTREAM_H_
#include <ostream>
#include <string>
#include <algorithm>
#include "BarHist.h"
#include "bar/BarAux.h"
#include "bar/BarCounter.h"
#include "bar/BarHist.h"

namespace ranger {
namespace bar {

void outputAll(const std::map<std::string, boost::icl::BarCounter>& counter,
        std::ostream& os);
void outputAll(const std::map<std::string, BarHist>& counter, std::ostream& os);
void outputAllWithHeader(const std::map<std::string, BarHist>& counter,
        const std::string& file_prefix, bool gzip=false);
void outputCounter(const boost::icl::BarCounter& counter,
        const std::string& chr, std::ostream& os);
void outputAllByChrWithHeader(std::map<std::string, BarHist>& counter,
        const std::string& file_prefix, const std::string& file_postfix, bool gzip=false);
void outputBar(const BarHist& hist, const std::string& chr, std::ostream& os);
void combineAndOutput(const BarHist& pos, const BarHist& neg,
        const std::string& chr, std::ostream& os);
void combineAndOutput(const BarHist& pos, const BarHist& neg,
        const std::string& chr, const std::string& file, bool gzip=false);
void combineAndOutputByChr(const std::map<std::string, BarHist>& pos,
        const std::map<std::string, BarHist>& neg, const std::string& chr,
        const std::string& file, bool gzip=false);

} /* namespace bar */
} /* namespace ranger */
#endif /* BGHISTSTREAM_H_ */
