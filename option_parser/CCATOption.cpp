/*
 * CCATOption.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#include "option_parser/CCATOption.h"
#include <string>
#include <iostream>
using namespace options::aux;
using namespace std;
using namespace boost;
using namespace boost::program_options;
namespace options {
int CCATOption::min_args = 2;
CCATOption::CCATOption() :
        version("1.16"), all("\nccat " + version + "\n\nUsage"), popt() {

    opt other("Other");
    opt input("Input");
    opt output("Output");
    opt qualities("Qualities");
    other.add_options()

    ("help,h", "show the usage")

    ("verbose", "show progress")

    ("version", "output the version number");

    input.add_options()

    ("data,d", po::value<string>(), "data file")

    ("control,c", po::value<string>(), "control file")

    ("format", po::value<string>(),
            "the format of the data file, can be one of : "
                    "bowtie, eland, sam, bam and bed");

    output.add_options()

    ("output,o", po::value<string>(), "the output location")

    ("report",

    "generate html reports")

    ("plot_region", po::value<uint32_t>()->default_value(6000),
            "the length of the snapshort regions in the report.")

    ("gene_annot_file", po::value<string>(),
            "the length of the snapshort regions in the report.");

    qualities.add_options()

    ("pval,p", po::value<double>()->default_value(0.0001), "p value cut-off")

    ("FDR,q", po::value<double>()->default_value(0.01), "FDR cut-off")

    ("ext_length,l", po::value<uint32_t>()->default_value(100),
            "read extension length");

    running_modes.add_options()

    ("thread,t", po::value<uint32_t>()->default_value(1),
            "number of running threads");

    popt.add("data", 1).add("control",1).add("output", 1);
    all.add(input).add(output).add(qualities).add(other);

}

CCATOption::~CCATOption() {

}

void CCATOption::parse(int _ac, char** _av) {
    hasEnoughArgs(_ac);
    po::store(
            po::command_line_parser(_ac, _av).options(all).positional(popt).run(),
            mVM);
    if (mVM.count("help")) {
        printHelp(all, cout);
        exit(0); //todo: should we move this up?
    }
    if (mVM.count("version")) {
        cout << "\n" << "1.16\n";
        exit(0); //todo: should we move this up?
    }
    verifyOptions();
}

#define VALUEOF(id) mVM[id].as<string>()

std::string CCATOption::printParsedOpts() {
    stringstream ss;
    ss << "\nprogram version:            " << version;
    ss << "\nData files:\n";
    ss << "\n File format:            " << mVM["format"].as<string>();
    ss << "\n Sample file:            " << mVM["data"].as<string>();
    ss << "\n Control file:           " << mVM["control"].as<string>();
    ss << "\nQualities:\n";
    ss << "\n P value cut off:        " << mVM["pval"].as<string>();
    ss << "\n FDR cut off:            " << mVM["FDR"].as<string>();
    ss << "\n Read extension length:  " << mVM["ext_length"].as<string>();
    ss << "\n Smoothing bandwidth:    " << mVM["bandwidth"].as<string>();
    ss << "\n Delta:                  " << mVM["delta"].as<string>() << endl;
    ss << "\n Pad region profile:     ";
    if (mVM.count("pad")) {
        ss << "Enabled";
    } else {
        ss << "Disabled\n";
    }
    ss << "Running modes:\n";
    ss << " Number of threads:      " << mVM["thread"].as<string>();
    ss << "\nOutput:\n";
    ss << " Regions:                "
            << mVM["out"].as<string>() + "_region.bed";
    ss << "\n Summits:                "
            << mVM["out"].as<string>() + "_summit.bed";
    ss << "\n Details of regions:     " << mVM["out"].as<string>() + "_details";
    ss << "\n HTML reports:           ";
    if (mVM.count("report")) {
        ss << "Enabled";
        ss << " Plot region length:     " << mVM["plot_region"].as<string>();
        ss << " Annotation file:        "
                << mVM["gene_annot_file"].as<string>();
    } else {
        ss << "Disabled(--report not specified)\n";
    }
    return ss.str();
}

void CCATOption::hasEnoughArgs(int argc) {
    if (argc < CCATOption::min_args) {
        printHelp(all, cout);
        throw std::logic_error("Not enough command options.");
    }
}

void CCATOption::verifyOptions() {
    require("data", mVM);
    require("output", mVM);
    require("format", mVM);

    file_r_good(mVM["data"].as<string>().c_str());
    file_w_good(mVM["output"].as<string>().c_str());

}

} /* namespace options */
