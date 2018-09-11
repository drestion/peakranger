/*
 * WigPEOption.cpp
 *
 *  Created on: Jun 25, 2012
 *      Author: xfeng
 */
#include "option_parser/OptionAux.h"
#include "option_parser/WigPEOption.h"
#include "common/boost_header.h"
#include "common/stl_header.h"


using namespace options::aux;
using namespace std;
using namespace boost;
using namespace boost::program_options;
namespace options {

int WigPEOption::min_args = 2;

WigPEOption::WigPEOption(const std::string& version) :
        version(version), all("\nWigPE " + version + "\n\nUsage"), popt() {
    opt other("Other");
    opt input("Input");
    opt output("Output");
    opt qualities("Qualities");
    other.add_options()

    ("help,h", "show the usage")

    ("verbose", "show progress when possible")

    ("version", "output the version number");

    input.add_options()

    ("data,d", po::value<string>(), "the data file");

    output.add_options()

    ("output,o", po::value<string>(), "the output location")

    ("split,s", "generate one wig file per chromosome")

    ("gzip,z", "compress the output")

    ("strand,x", "generate one wig file per strand");

    qualities.add_options()

    ("ext_length,l", po::value<uint32_t>()->default_value(0),
            "read extension length");

    popt.add("data", 1).add("output", 1);
    all.add(input).add(output).add(qualities).add(other);

}

WigPEOption::~WigPEOption() {

}

void WigPEOption::parse(int _ac, char** _av) {
    hasEnoughArgs(_ac);
    po::store(
            po::command_line_parser(_ac, _av).options(all).positional(popt).run(),
            mVM);
    po::notify(mVM);
    if (mVM.count("help")) {
        printHelp(all, cout);
        exit(0); //todo: should we move this up?
    }
    if (mVM.count("version")) {
        printVersion(cout, version);
    }
    verifyOptions();
}

string WigPEOption::printParsedOpts() {
    stringstream os;

    os << "\n" << "Wig version:              " << version;
    os << "\n" << "Data files:";
    os << "\n" << " Sample file:             " << mVM["data"].as<string>();
    os << "\n" << "Qualities:";
    os << "\n" << " Read extension length:   "
            << mVM["ext_length"].as<uint32_t>();
    os << "\n" << "Output:";
    os << "\n" << " Gzip results:            ";
    if (mVM.count("gzip")) {
        os << "Yes";
    } else {
        os << "No";
    }
    os << "\n" << " One wig per strand:      ";
    if (mVM.count("strand")) {
        os << "Yes";
    } else {
        os << "No";
    }
    if (mVM.count("split")) {
        os << "\n" << " Splitting results:       Yes";
        os << "\n" << " Result file:             "
                << mVM["output"].as<string>();
        os << "\n" << "                          "
                << "and other splitted files in this directory";
    } else {
        os << "\n" << " Splitting results:       No";
        os << "\n" << " Result file:             "
                << mVM["output"].as<string>();
    }
    return os.str();
}

void WigPEOption::hasEnoughArgs(int argc) {
    if (argc < WigPEOption::min_args) {
        printHelp(all, cout);
        throw std::logic_error("Not enough command options.");
    }
}

void WigPEOption::verifyOptions() {
    require("data", mVM);
    require("output", mVM);
    file_r_good(mVM["data"].as<string>().c_str());
    file_w_good(mVM["output"].as<string>().c_str());
}

} /* namespace options */
