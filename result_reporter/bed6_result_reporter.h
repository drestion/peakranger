/*
 * bed6_result_reporter.h
 *
 *  Created on: Jun 8, 2011
 *      Author: xin
 */

#ifndef BED6_RESULT_REPORTER_H_
#define BED6_RESULT_REPORTER_H_
#include <ostream>
#include "result_reporter.h"
#include <stdint.h>
#include <vector>
#include <map>

typedef std::map<std::string, std::vector<called_peak> > enriched_regions;

class bed6_result_reporter: public result_reporter {
public:
    bed6_result_reporter();
    ~bed6_result_reporter();

    void report_fdr_summit(enriched_regions& result, std::ostream& om);
    void report_pval_region(enriched_regions& result, std::ostream& om);
    void report_pval_summit(enriched_regions& result, std::ostream& om);
    void report_fdr_region(enriched_regions& regions, std::ostream& om);
};

#endif /* BED6_RESULT_REPORTER_H_ */
