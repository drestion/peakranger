/*
 * AppAux.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef APPAUX_H_
#define APPAUX_H_
#include <string>

#include "common/boost_header.h"
#include "parser/readsParser.h"
#include "utils/exceptions.h"
#include "parser/bowtieParser.h"
#include "parser/samParser.h"
#include "parser/bamParser.h"
#include "parser/bedParser.h"
#include "option_parser/cmd_option_parser.h"
#include "option_parser/cmd_option_parser.h"
#include "region_detector/region_detector.h"
#include "region_detector/calledpeak.h"
#include "region_detector/fdr_based_thresholder.h"
#include "region_detector/ccat.h"
#include "result_reporter/result_reporter.h"
#include "result_reporter/bed6_result_reporter.h"
#include "wiggle/wiggle_reporter.h"
#include "wiggle/JTwigglefile.h"

#include "utils/logger.h"
#include "utils/timer.h"
#include "utils/stringutil.h"
#include "utils/util_print.h"
#include "utils/Tracer.h"
#include "short_reads/readstools.h"

#include "ggplay/chipseqhtmlreporter.h"
#include "region_detector/ccat_main.h"
#include "app/AppAux.h"
namespace app {
namespace aux {

void CalculateFDR(cmd_option_parser & option,
        boost::shared_ptr<region_detector> & detector);
void OutputResults(cmd_option_parser & option,
        boost::shared_ptr<region_detector> & detector,
        const size_t & fdr_passed, const size_t & fdr_failed);

void getParser(cmd_option_parser & option,
        boost::shared_ptr<readsParser> & parser);

void parseReads(cmd_option_parser& option, const std::string& readsfile,
        boost::shared_ptr<readsParser>& parser, Reads& treads);

void filterByFDR(const std::map<std::string, std::vector<called_peak> >& toFilter,
        std::map<std::string, std::vector<called_peak> >& results, double fdr);

}
} /* namespace app */
#endif /* APPAUX_H_ */
