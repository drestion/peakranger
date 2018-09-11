/*
 * peakranger_cmd_option_parser
 *
 *  Created on: July 25, 2011
 *      Author: xin
 */

#include "peakranger_cmd_option_parser.h"
#include "utils/stringutil.h"
#include "utils/util_print.h"
#include "option_parser/OptionAux.h"
#include <stdlib.h>
#include <stdint.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "utils/exceptions.h"
using namespace std;
using namespace boost;
using namespace utils;
using namespace boost::program_options;
using namespace options::aux;
using namespace options;
peakranger_cmd_option_parser::peakranger_cmd_option_parser(int argc,
        char** argv, const std::string& version) :
        cmd_option_parser(argc, argv), version(version), all(
                "\nranger " + version + "\n\nUsage"), popt() {
    _binlength = 10000;
    opt other("Other");
    opt input("Input");
    opt output("Output");
    opt qualities("Qualities");
    opt running_modes("Running modes");

    other.add_options()

    ("help,h", "show the usage")

    ("verbose", "show progress")

    ("version", "output the version number");

    input.add_options()

    ("data,d", po::value<string>(&_treat_dir), "data file")

    ("control,c", po::value<string>(&_control_dir), "control file")

    ("format", po::value<string>(&_format),
            "the format of the data file, can be one of : "
                    "bowtie, sam, bam and bed");

//	("config", po::value<string>(&_config_file), "configuration file");

    output.add_options()

    ("output,o", po::value<string>(&_output_dir), "the output location")

    ("report",

    "generate html reports")

    ("plot_region",
            po::value<uint32_t>(&_html_region_length)->default_value(6000),
            "the length of the snapshort regions in the report")

    ("gene_annot_file", po::value<string>(&_gene_anno_file),
            "the gene annotation file");

    qualities.add_options()

    ("pval,p", po::value<double>(&_p_cut_off)->default_value(double(0.0001L)),
            "p value cut-off")

    ("FDR,q", po::value<double>(&_fdr_cut_off)->default_value(double(0.01L)),
            "FDR cut-off")

    ("ext_length,l", po::value<uint32_t>(&_ext_length)->default_value(100),
            "read extension length")

    ("delta,r", po::value<double>(&_delta)->default_value(double(1)),
            "sensitivity of the summit detector")

    ("bandwidth,b", po::value<uint32_t>(&_bandwidth)->default_value(99),
            "smoothing bandwidth")

    ("pad", "pad read coverage profile to avoid false positive summmits");

    running_modes.add_options()("thread,t",
            po::value<uint32_t>(&_no_of_thread)->default_value(maxThreads - 1),
            "number of worker threads");
    _config_file = "";

    popt.add("data", 1).add("control", 1).add("output", 1);
    all.add(input).add(output).add(qualities).add(running_modes).add(other);
}

peakranger_cmd_option_parser::~peakranger_cmd_option_parser() {
}
void peakranger_cmd_option_parser::parse() {

    if (_ac < 2) {
        options::aux::printHelp(all, cout);

        throw std::logic_error(
                "\nNot enough command options.\n\nProvided args:\n"
                        + options::aux::printRawOpts(_ac, _av));
    }

    po::store(
            po::command_line_parser(_ac, _av).options(all).positional(popt).run(),
            vm);
    po::notify(vm);

    setFormat(to_lower_copy(trim_copy(getFormat())));

    if (vm.count("help")) {
        setHelpRequested(true);
        options::aux::printHelp(all, cout);
        exit(0);
    }
    if (vm.count("verbose")) {
        setVerboseRequested(true);
    }

    if (vm.count("version")) {
        printVersion(cout, version);

    }

    if (vm.count("report")) {
        this->setNeedHtml(true);
    } else {
        this->setNeedHtml(false);
    }
    if (vm.count("split")) {
        this->setSplit(true);
    } else {
        this->setSplit(false);
    }

    if (vm.count("gene_annot_file")) {
        this->setNeedHtml(true);
        file_r_good(_gene_anno_file.c_str());
    }
    require("data", vm);
    require("control", vm);
    require("output", vm);
    require("format", vm);
    if (vm.count("report")) {
        require("gene_annot_file", vm);
    }
    file_r_good(_treat_dir.c_str());
    file_r_good(_control_dir.c_str());
    string dir, file, file_ext;

    stringutil::get_dir_file(_treat_dir, dir, file, file_ext);

    setTreat_file(_treat_dir);
    setTreat_dir(dir);
    setTreatfilename(file);

    stringutil::get_dir_file(_control_dir, dir, file, file_ext);

    setControl_file(_control_dir);
    setControl_dir(dir);
    setControlfilename(file);

    file_w_good(_output_dir.c_str());
    stringutil::get_dir_file(_output_dir, dir, file, file_ext);
    //todo: linux only
    setTreat_wig_file(dir + "/" + _treatfilename + ".wig");
    setControl_wig_file(dir + "/" + _controlfilename + ".wig");
    setOutput_file(_output_dir);
    setOutput_dir(dir);

    setReportName(dir + "/reports");

    if (vm.count("pad")) {
        setPad(true);
    } else {
        setPad(false);
    }

    verify();

}

void peakranger_cmd_option_parser::print_option(ostream & os) {
    os << ("program version:          ") << version << endl;
    os << ("Data files:\n");
    os << (" File format:             ") << getFormat() << endl;
    os << (" Sample file:             ") << getTreat_file() << endl;

    os << (" Control file:            ") << getControl_file() << endl;

    os << ("Qualities:\n");
    os << (" P value cut off:         ") << getCut_off() << endl;
    os << (" FDR cut off:             ") << getFdrCutOff() << endl;
    os << (" Read extension length:   ") << getExt_length() << endl;
    os << (" Smoothing bandwidth:     ") << getBandwidth() << endl;
    os << (" Delta:                   ") << getDelta() << endl;
    os << (" Pad region profile:      ");
    if (getPad()) {
        os << "Enabled" << endl;
    } else {
        os << "Disabled\n";
    }
    os << ("Running modes:\n");

    os << (" Number of threads:       ") << getNo_of_thread() << endl;
    os << ("Output:\n");

    os << (" Regions:                 ") << getOutput_file() + "_region.bed"
            << endl;
    os << (" Summits:                 ") << getOutput_file() + "_summit.bed"
            << endl;
    os << (" Details of regions:      ") << getOutput_file() + "_details"
            << endl;
    os << (" HTML reports:            ");
    if (needHtml()) {
        os << "Enabled" << endl;
        os << (" Plot region length:      ") << getHtmlRegionLength() << endl;
        os << (" Annotation file:         ") << getGeneAnnoFile() << endl;
    } else {
        os << "Disabled(--report not specified)\n";
    }

}

bool peakranger_cmd_option_parser::isSplit() const {
    return _split;
}

void peakranger_cmd_option_parser::setSplit(bool _split) {
    this->_split = _split;
}

void peakranger_cmd_option_parser::print_option_file(ostream & os) const {
    os << ("#program version:           ") << version << endl;
    os << ("#Data files:\n");
    os << ("# File format:             ") << getFormat() << endl;
    os << ("# Sample file:             ") << getTreat_file() << endl;
    os << ("# Control file:            ") << getControl_file() << endl;
    os << ("#Qualities:\n");
    os << ("# P value cut off:         ") << getCut_off() << endl;
    os << ("# FDR cut off:             ") << getFdrCutOff() << endl;
    os << ("# Read extension length:   ") << getExt_length() << endl;
    os << ("# Smoothing bandwidth:     ") << getBandwidth() << endl;
    os << ("# Delta:                   ") << getDelta() << endl;
    os << ("# Pad region profile:      ");
    if (getPad()) {
        os << "Enabled" << endl;
    } else {
        os << "Disabled\n";
    }
    os << ("#Running modes:\n");
    os << ("# Number of threads:       ") << getNo_of_thread() << endl;
    os << ("#Output:\n");

    os << ("# Regions:                 ") << getOutput_file() + "_region.bed"
            << endl;
    os << ("# Summits:                 ") << getOutput_file() + "_summit.bed"
            << endl;
    os << ("# Details of regions:      ") << getOutput_file() + "_details"
            << endl;
    os << ("# HTML reports:            ");
    if (needHtml()) {
        os << "Enabled" << endl;
        os << ("# Plot region length:      ") << getHtmlRegionLength() << endl;
        os << ("# Annotation file:         ") << getGeneAnnoFile() << endl;
    } else {
        os << "Disabled(--report not specified)\n";
    }

}

bool peakranger_cmd_option_parser::needHtml() const {
    return _html;
}

void peakranger_cmd_option_parser::setNeedHtml(bool _html) {
    this->_html = _html;
}

void peakranger_cmd_option_parser::verify() {
    std::ostringstream oss;
    oss << "Specified ";
    if (!is_in_range(_delta, 0.0, 1.000001)) {
        oss << "delta value " << _delta << " is not in the range : (0, 1)"
                << endl;
        throw not_in_range(oss.str().c_str());
    }

    if (!is_in_range<double>(_p_cut_off, 0.0, 1.0)) {
        oss << "p-value cut off " << _p_cut_off
                << " is not in the range : (0, 1)" << endl;
        throw not_in_range(oss.str().c_str());
    }

    if (!is_in_range<uint32_t>(_html_region_length, 1, 1000000)) {
        oss << "region plot length " << _html_region_length
                << " is not in the range : (1, 1000000)" << endl;
        throw not_in_range(oss.str().c_str());
    }

    if (_format != "bed" && _format != "bowtie" && _format != "eland"
            && _format != "sam" && _format != "bam") {
        oss << "file format " << _format
                << " is not a valid entry. Please choose from : bowtie, sam, bam and bed."
                << endl;
        throw not_in_range(oss.str().c_str());
    }


    if (maxThreads < getNo_of_thread()) {
        oss << "number of threads " << getNo_of_thread()
                << " is not valid. The current system allows up to "
                << maxThreads << " threads.\n";
        throw not_in_range(oss.str().c_str());
    }

}

