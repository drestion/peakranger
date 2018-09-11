/*
 * bed6_result_reporter.cpp
 *
 *  Created on: Jun 8, 2011
 *      Author: xin
 */

#include "bed6_result_reporter.h"
#include <string>
#include <map>
#include <vector>
#include <stdint.h>
#include <utility>
#include <sstream>
using namespace std;

typedef map<string, vector<called_peak> > enriched_regions;

bed6_result_reporter::bed6_result_reporter() {

}

bed6_result_reporter::~bed6_result_reporter() {

}

void bed6_result_reporter::report_fdr_summit(enriched_regions & regions,
                                             std::ostream & om)
                                             {
    enriched_regions::iterator it = regions.begin();
    string chr;
    print_cite(om);
    om << "#summit_chr\tsummit_start\tsummit_end\tsummit_ID\tsummit_FDR\tsummit_strand\n";
    uint64_t cnt = 0;

    for (; it != regions.end(); it++) {

        chr = (it->first);

        vector<called_peak>::iterator rit = it->second.begin();
        for (; rit != it->second.end(); rit++) {
            std::ostringstream s;
            s << "PeakRanger_" << cnt++ << "_region_" << rit->first << "_"
            << rit->second << "_pval_" << rit->p << "_FDR_"
            << rit->q;
            vector<uint32_t>::iterator sit = rit->summits.begin();
            for (; sit != rit->summits.end(); sit++) {
                om << chr << "\t" << *sit << "\t" << (*sit) + 1 << "\t"
                << s.str() << "\t" << rit->q << "\t+\n";
            }
        }
    }
}

void bed6_result_reporter::report_pval_region(enriched_regions & regions,
                                              std::ostream & om)
                                              {
    enriched_regions::iterator it = regions.begin();
    string chr;
    print_cite(om);
    om << "#region_chr\tregion_start\tregion_end\tregion_ID\tregion_pvalue\tregion_strand\n";
    uint64_t cnt = 0;
    for (; it != regions.end(); it++) {

        chr = (it->first);

        vector<called_peak>::iterator rit = it->second.begin();
        for (; rit != it->second.end(); rit++) {
            std::ostringstream s;
            s << "PeakRanger_" << cnt++ << "_pval_" << rit->p << "_FDR_"
            << rit->q;
            om << chr << "\t" << rit->first << "\t" << rit->second << "\t"
            << s.str() << "\t" << rit->p << "\t+\n";
        }
    }
}

void bed6_result_reporter::report_pval_summit(enriched_regions & regions,
                                              std::ostream & om)
                                              {
    cout <<"Not implemented yet.\n";
}

void bed6_result_reporter::report_fdr_region(enriched_regions& regions,
                                             std::ostream& om) {
    enriched_regions::iterator it = regions.begin();
    string chr;
    print_cite(om);
    om << "#region_chr\tregion_start\tregion_end\tregion_ID\tregion_FDR\tregion_strand\n";
    uint64_t cnt = 0;

    for (; it != regions.end(); it++) {

        chr = (it->first);

        vector<called_peak>::iterator rit = it->second.begin();
        for (; rit != it->second.end(); rit++) {
            std::ostringstream s;
            s << "PeakRanger_" << cnt++ << "_pval_" << rit->p << "_FDR_"
            << rit->q;
            om << chr << "\t" << rit->first << "\t" << rit->second << "\t"
            << s.str() << "\t" << rit->q << "\t+\n";
        }
    }
}

