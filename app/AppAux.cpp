/*
 * AppAux.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#include <stdexcept>
#include <fstream>
#include "common/boost_header.h"
#include "app/AppAux.h"
using namespace std;
using namespace utils;
using namespace boost;
typedef map<string, vector<called_peak> > enriched_regions;
typedef vector<called_peak>::iterator ritrr;
typedef enriched_regions::iterator pritrr;

namespace app {
namespace aux {

void CalculateFDR(cmd_option_parser & option,
        boost::shared_ptr<region_detector> & detector) {
//    if (option.getVerboseRequested()) {
//        cout << "\n Calculating FDR...\n" << "\n";
//    }
//    string chr;
//    pritrr it = detector->_resultRegions.begin();
//    for (; it != detector->_resultRegions.end(); it++) {
//        chr = it->first;
//        size_t d_pk_cnt = it->second.size();
//        if (d_pk_cnt < 1) {
//            continue;
//        }
//        std::sort(detector->_resultRegions[chr].begin(),
//                detector->_resultRegions[chr].end(), sorter_by_pval);
//        size_t j = 0;
//        size_t _rk = 0;
//        vector<called_peak> & _pk = detector->_resultRegions[chr];
//        for (; j < _pk.size(); j++) {
//            if (_pk[j].p < 0L) {
//                _pk[j].q = 0;
//            } else {
//                _rk++;
//                _pk[j].q = _pk[j].p * _pk.size() / _rk;
//            }
//        }
//
//    }

}

void OutputResults(cmd_option_parser & option,
        boost::shared_ptr<region_detector> & detector, size_t & fdr_passed,
        size_t & fdr_failed) {
    string ga;
    ga = (option.getOutput_file() + "_region.bed");
    ofstream of(ga.c_str());
    if (!(of.is_open())) {
        throw FileNotGood(ga.c_str());
    }
    ga = (option.getOutput_file() + "_summit.bed");
    ofstream of_smt(ga.c_str());
    if (!(of_smt.is_open())) {
        throw FileNotGood(ga.c_str());
    }
    ga = (option.getOutput_file() + "_details");
    ofstream of_raw(ga.c_str());
    if (!(of_raw.is_open())) {
        throw FileNotGood(ga.c_str());
    }
    utilprint::citation ct;
    ct.print_msg(of_raw);
    ct.print_msg(of_smt);
    ct.print_msg(of);
    logDate(of_raw);
    logDate(of_smt);
    logDate(of);
    option.print_option_file(of_raw);
    option.print_option_file(of_smt);
    option.print_option_file(of);
    of_smt
            << "\n#summit_chr\tsummit_start\tsummit_end\tsummit_ID\tsummit_FDR\tsummit_strand\n";
    of
            << "\n#region_chr\tregion_start\tregion_end\tregion_ID\tregion_pvalue\tregion_strand\n";
    of_raw
            << "\n#region_chr\tregion_start\tregion_end\tregion_ID\tregion_summits\tregion_pvalue\tregion_fdr\tregion_strand\tregion_treads\tregion_creads\n";

    pritrr it = detector->_resultRegions.begin();
    for (; it != detector->_resultRegions.end(); it++) {

        foreach (called_peak pk , it->second) {
            of_raw << it->first << "\t" << pk.first << "\t" << pk.second
                    << "\tranger";
            if (pk.q <= option.getFdrCutOff()) {
                of_raw << "_fdrPassed_" << fdr_passed;
            } else {
                of_raw << "_fdrFailed_" << fdr_failed;
            }
            of_raw << "_pval_" << pk.p << "_fdr_" << pk.q;
            of_raw << "\t";
            vector<uint32_t>::iterator sit = pk.summits.begin();
            for (; sit != pk.summits.end(); sit++) {
                of_raw << *sit << "_";
                of_smt << it->first << "\t" << *sit << "\t" << (*sit) + 1
                        << "\tranger_region_" << pk.first << "_" << pk.second
                        << "_pval_" << pk.p;
                if (pk.q <= option.getFdrCutOff()) {
                    of_smt << "_fdrPassed_" << fdr_passed;
                } else {
                    of_smt << "_fdrFailed_" << fdr_failed;
                }

                of_smt << "\t" << pk.q << "\t+\n";
            }

            of_raw << "\t" << pk.p << "\t" << pk.q << "\t+\t" << pk.treads
                    << "\t" << pk.creads << "\n";

            of << it->first << "\t" << pk.first << "\t" << pk.second
                    << "\tranger";
            if (pk.q <= option.getFdrCutOff()) {
                of << "_fdrPassed_" << fdr_passed++;
            } else {
                of << "_fdrFailed_" << fdr_failed++;
            }
            of << "_pval_" << pk.p << "_fdr_" << pk.q << "\t" << pk.q
                    << "\t+\n";
        }
    }
}

void getParser(cmd_option_parser & option,
        boost::shared_ptr<readsParser> & parser) {
    if (option.getFormat() == cmd_option_parser::format_bowtie) {
        parser = boost::make_shared<bowtieParser>();
    } else if (option.getFormat() == cmd_option_parser::format_sam) {
        parser = boost::make_shared<samParser>();
    } else if (option.getFormat() == cmd_option_parser::format_bed) {
        parser = boost::make_shared<bedParser>();
    } else if (option.getFormat() == cmd_option_parser::format_bam) {
        parser = boost::make_shared<bamParser>();
    } else {
        string str("The specified format ");
        str += option.getFormat();
        str += " has not been implemented yet.\n";
        throw not_in_range(str.c_str());
    }

}

void parseReads(cmd_option_parser& option, string readsfile,
        boost::shared_ptr<readsParser>& parser, Reads& treads) {

    if (option.getChrtableSpecified()) {
        vector<string> chrs_to_parse = option.getChrs_to_parse();
        parser->parse(treads, readsfile, chrs_to_parse);
    } else {
        parser->parse(treads, readsfile);
    }

}

void filterByFDR(const map<string, vector<called_peak> >& toFilter,
        map<string, vector<called_peak> >& passFDR, double fdr) {

    map<string, vector<called_peak> >::const_iterator it;
    it = toFilter.begin();
    for (; it != toFilter.end(); it++) {
        foreach (called_peak pk , it->second) {
            if (pk.q <= fdr) {
                passFDR[it->first].push_back(pk);
            }
        }
    }

}
}
} /* namespace app */
