/*
 * LCOption.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: xfeng
 */

#include "option_parser/LCOption.h"
#include "option_parser/OptionAux.h"
#include "common/boost_header.h"
#include "common/stl_header.h"

using namespace options::aux;
using namespace std;
using namespace boost;
using namespace boost::program_options;
namespace options {

LCOption::~LCOption() {

}

int LCOption::min_args = 2;

LCOption::LCOption(const std::string& version) :
        version(version), all("\nlc " + version + "\n\nUsage"), popt() {
    opt other("Other");
    opt input("Input");

    other.add_options()

    ("help,h", "show the usage")

    ("verbose", "show progress when possible")

    ("version", "output the version number");

    input.add_options()

    ("data,d", po::value<string>(), "the data file");

    popt.add("data", 1);
    all.add(input).add(other);

}

void LCOption::parse(int _ac, char** _av) {
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

string LCOption::printParsedOpts() {
    stringstream os;

    os << "\n" << "lc version:               " << version;
    os << "\n" << " Sample file:             " << mVM["data"].as<string>();

    return os.str();
}

void LCOption::hasEnoughArgs(int argc) {
    if (argc < LCOption::min_args) {
        printHelp(all, cout);
        throw std::logic_error("Not enough command options.");
    }
}

void LCOption::verifyOptions() {
    require("data", mVM);
    file_r_good(mVM["data"].as<string>().c_str());
}

} /* namespace options */
