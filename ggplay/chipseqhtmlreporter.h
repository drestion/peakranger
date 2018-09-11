/*
 * chipseqhtmlreporter.h
 *
 *  Created on: Feb 15, 2012
 *      Author: xfeng
 */

#ifndef CHIPSEQHTMLREPORTER_H_
#define CHIPSEQHTMLREPORTER_H_

#include "wiggle/wigbuilder.h"
#include "region_detector/calledpeak.h"
#include "option_parser/cmd_option_parser.h"
#include "tab_file/TabGene.h"
#include "tab_file/Region.h"
#include "tab_file/NearbyGeneFinder.h"
#include <iostream>
#include <string>
#include <stdint.h>
#include <map>
#include "short_reads/reads.h"

//todo: refactor
class chipseq_html_reporter {
public:
    chipseq_html_reporter();

    void generate_report(Reads& treads, Reads& creads,
            std::map<std::string, std::vector<called_peak> >& peaks,
            cmd_option_parser& option);

    void setRegionLength(uint32_t _d);
    void addCitation(const std::string& citation);
private:
    void printCitation(std::ostream& os);
    void print_chr(std::string& chr,
            std::map<std::string, std::vector<called_peak> >& peaks,
            std::ostream& os);
    void print_peak(called_peak& pk, std::string chr, std::ostream& os);
    void print_head(std::ostream& os,uint32_t reportLength);
    void print_end(std::ostream& os);
    void print_css(cmd_option_parser& option);

    void prepare_dir(cmd_option_parser& option, const std::string& reportname);
    void run_img_scripts(cmd_option_parser& option,
            std::map<std::string, std::vector<called_peak> > & peaks);
    void print_html_table(
            std::map<std::string, std::vector<called_peak> > & peaks,
            cmd_option_parser & option, std::ostream& os);
    void print_pk_img_wig(wigs& wt, wigs& wc, called_peak& pk,
            std::ostream& os);
    void print_pk_img_para(std::ostream& os, called_peak& pk, std::string& rf);
    void print_peak_img_script(called_peak & pk, std::string rf, wigs& wt,
            wigs& wc, std::ostream& os, std::string chr);
    void prepare_wigs(Reads& treads, Reads& creads, std::string& chr,
            cmd_option_parser& option, wigs& wt, wigs& wc);
    void print_img_script(Reads& treads, Reads& creads,
            std::map<std::string, std::vector<called_peak> > & peaks,
            cmd_option_parser& option);
    void print_r_head(std::ostream& os);
    void print_r_tc_drawing(std::ostream& os);
    void run_pk_img_script(std::string& rf);

    void format_name();
    void print_gene_annot(const called_peak& pk, const std::string& chr,
            std::ostream& os);
    void prepare_gene_annot();
    void print_axis_to_end(std::ostream& os);
    void prepare_report_name(cmd_option_parser & option);
    void filter_in_range_anno(std::vector<tab_file::TabGene> & genes_to_filter,
            std::vector<tab_file::TabGene> & genes_to_print);
    void get_genes_to_filter(const std::string & chr, tab_file::Region & pkrg,
            std::vector<tab_file::TabGene> & genes_to_filter);
    wig_builder _wig;
    uint32_t _d;
    uint32_t _binlength;
    std::string _report_name;
    std::string _gene_anno_file;
    std::map<std::string, std::vector<tab_file::TabGene> > _genes;
    tab_file::NearbyGeneFinder _nbgf;
    std::string _cite;
};

#endif /* CHIPSEQHTMLREPORTER_H_ */
