/*
 * NROption.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#include "option_parser/NROption.h"
#include <string>
#include <exception>
using namespace options::aux;
using namespace std;
using namespace boost;
using namespace boost::program_options;

namespace options {
int NROption::min_args = 2;
NROption::NROption(const string& version) :
		version(version), all("\nnr " + version + "\n\nUsage"), popt() {
	opt other("Other");
	opt input("Input");
	opt output("Output");
	opt qualities("Qualities");
	other.add_options()

	("help,h", "show the usage")

	("verbose", "show progress when possible")

	("version", "output the version number");

	input.add_options()

	("data,d", po::value<string>(), "data file")

	("control,c", po::value<string>(), "control file")

	("format", po::value<string>(),
			"the format of the data file, can be one of : "
					"bowtie, sam, bam and bed");
	qualities.add_options()

	("ext_length,l", po::value<uint32_t>()->default_value(200),
			"read extension length");

	popt.add("data", 1).add("control", 1);
	all.add(input).add(qualities).add(other);

}

NROption::~NROption() {

}

void NROption::parse(int _ac, char** _av) {
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
		cout << "\n" << version << "\n\n";
		exit(0); //todo: should we move this up?
	}
	verifyOptions();
}

std::string NROption::printParsedOpts() {
	stringstream ss;

	ss << "\n" << "program version:          " << version << "\n";
	ss << "\n" << "Data files:";
	ss << "\n" << " file format:             " << mVM["format"].as<string>();
	ss << "\n" << " Sample file:             " << mVM["data"].as<string>();
	ss << "\n" << " Control file:            " << mVM["control"].as<string>();
	ss << "\n" << "Qualities:            ";
	ss << "\n" << " Read extension:          "
			<< mVM["ext_length"].as<uint32_t>();
	ss << "\n";
	return ss.str();
}

void NROption::hasEnoughArgs(int argc) {
	if (argc < NROption::min_args) {
		printHelp(all, cout);
		throw std::logic_error("Not enough command options.");
	}
}

void NROption::verifyOptions() {
	require("data", mVM);
	require("control", mVM);
	require("format", mVM);
	string _format = mVM["format"].as<string>();
	if (_format != "bed" && _format != "bowtie" && _format != "sam"
			&& _format != "bam") {
		stringstream oss;
		oss << "file format " << _format
				<< " is not a valid entry. Please choose from : bowtie, sam, bam and bed."
				<< endl;
		throw std::logic_error(oss.str().c_str());
	}
	file_r_good(mVM["data"].as<string>().c_str());
}

} /* namespace options */
