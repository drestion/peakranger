/*
 * region_detector.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: xin
 */

#include "region_detector.h"
#include "utils/exceptions.h"
#include "utils/assert_helpers.h"
using namespace std;

//
//typedef map<string, vector<called_peak> > enriched_regions;
//typedef vector<called_peak>::iterator ritrr;
//typedef enriched_regions::iterator pritrr;

void region_detector::export_results(result_reporter & reporter,
                                     ostream & om) {
    cout <<"gagag*****************export"<<endl;
    reporter.report_fdr_region(_resultRegions,
                    om);
}
