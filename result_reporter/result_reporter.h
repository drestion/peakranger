/*
 * result_reporter.h
 *
 *  Created on: Jun 7, 2011
 *      Author: xin
 */

#ifndef RESULT_REPORTER_H_
#define RESULT_REPORTER_H_
#include <ostream>
#include <vector>
#include <utility>
#include <map>
#include <stdint.h>
#include <iostream>
#include <utility>

#include "region_detector/calledpeak.h"



typedef std::map<std::string, std::vector<called_peak> > enriched_regions;

class result_reporter {
public:
    result_reporter();
    virtual ~result_reporter();
    virtual void report_fdr_region(enriched_regions& result,
            std::ostream& om)=0;
    virtual void report_fdr_summit(enriched_regions& result,
            std::ostream& om)=0;
    virtual void report_pval_region(enriched_regions& result,
            std::ostream& om)=0;
    virtual void report_pval_summit(enriched_regions& result,
            std::ostream& om)=0;
    void print_cite(std::ostream& om);
};

#endif /* RESULT_REPORTER_H_ */
