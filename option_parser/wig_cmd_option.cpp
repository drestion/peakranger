/*
 * boostCMDOption.cpp
 *
 *  Created on: May 28, 2011
 *      Author: xin
 */

#include "wig_cmd_option.h"
#include "utils/stringutil.h"
#include "option_parser/OptionAux.h"
#include <stdlib.h>
#include <stdint.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "option_parser/OptionAux.h"
#include "utils/exceptions.h"
using namespace std;
using namespace boost;
using namespace utils;
using namespace boost::program_options;
namespace po = boost::program_options;
using namespace options::aux;

wig_cmd_option::~wig_cmd_option() {
}
wig_cmd_option::wig_cmd_option(int argc, char** argv,
		const std::string& version) :
		cmd_option_parser(argc, argv), version(version), all("\nwig " + version + "\n\nUsage"), _split(
				false), _gz(false), _stranded(false) {
	typedef options_description opt;
	typedef positional_options_description p_opt;

	opt other("Other");
	opt input("Input");
	opt output("Output");
	opt qualities("Qualities");

	other.add_options()

	("help,h", "show the usage")

	("verbose", "show progress")

	("version", "output the version number");

	input.add_options()

	("data,d", po::value<string>(&_treat_dir), "data file")

	("format", po::value<string>(&_format),
			"the format of the data file, can be one of : "
					"bowtie, sam, bam and bed");

	output.add_options()

	("output,o", po::value<string>(&_output_dir), "the output location")

	("split,s", "generate one wig file per chromosome")

	("gzip,z", "compress the output")

	("strand,x", "generate one wig file per strand");

	qualities.add_options()

	("ext_length,l", po::value<uint32_t>(&_ext_length)->default_value(200),
			"read extension length");

	popt.add("data", 1).add("output", 1);
	all.add(input).add(output).add(qualities).add(other);
}
void wig_cmd_option::parse() {
	if (_ac < 2) {
		options::aux::printHelp(all, cout);
		cout << "\nNot enough args\n";
		cout << "\nProvided args:\n" << printRawOpts(_ac, _av) << "\n";

		exit(0);
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
	if (vm.count("split")) {
		this->setSplit(true);
	} else {
		this->setSplit(false);
	}
	if (vm.count("strand")) {
		this->setStranded(true);
	} else {
		this->setStranded(false);
	}
	if (vm.count("gzip")) {
		this->setGz(true);
	} else {
		this->setGz(false);
	}
	require("data", vm);
	require("output", vm);
	require("format", vm);

	file_r_good(_treat_dir.c_str());
	string dir, file, file_ext;

	stringutil::get_dir_file(_treat_dir, dir, file, file_ext);

	setTreat_file(_treat_dir);
	setTreat_dir(dir);
	setTreatfilename(file);

	file_w_good(_output_dir.c_str());
	stringutil::get_dir_file(_output_dir, dir, file, file_ext);
	//todo: linux only
	if (file_ext != "wig") {
		_output_dir += ".wig";
	}
	setOutput_file(_output_dir);
	setOutput_dir(dir);

	verify();

}

void wig_cmd_option::print_option_file(ostream & os) const {
	os << "\n" << "program version:          " << version << "\n";
	os << ("Data files:\n");
	os << (" File format:             ") << getFormat() << endl;
	os << (" Sample file:             ") << getTreat_file() << endl;

	os << ("Qualities:\n");

	os << (" Read extension length:   ") << getExt_length() << endl;

	os << ("Output:\n");
	if (isGz()) {
		os << (" Gzip results:            Yes") << endl;
	} else {
		os << (" Gzip results:            No") << endl;
	}
	if (isStranded()) {
		os << (" One wig per strand:      Yes") << endl;
	} else {
		os << (" One wig per strand:      No") << endl;
	}
	if (isSplit()) {
		os << (" Splitting results:       Yes") << endl;
		os << (" Result file:             ") << getOutput_file() << endl;
		os
				<< ("                          and other splitted files in this directory")
				<< endl;
	} else {
		os << (" Splitting results:       No") << endl;
		os << (" Result file:             ") << getOutput_file() << endl;
	}
}
void wig_cmd_option::print_option(ostream & os) {
	os << "\n" << "program version:          " << version << "\n";
	os << ("Data files:\n");
	os << (" File format:             ") << getFormat() << endl;
	os << (" Sample file:             ") << getTreat_file() << endl;

	os << ("Qualities:\n");

	os << (" Read extension length:   ") << getExt_length() << endl;

	os << ("Output:\n");
	if (isGz()) {
		os << (" Gzip results:            Yes") << endl;
	} else {
		os << (" Gzip results:            No") << endl;
	}
	if (isStranded()) {
		os << (" One wig per strand:      Yes") << endl;
	} else {
		os << (" One wig per strand:      No") << endl;
	}
	if (isSplit()) {
		os << (" Splitting results:       Yes") << endl;
		os << (" Result file:             ") << getOutput_file() << endl;
		os
				<< ("                          and other splitted files in this directory")
				<< endl;
	} else {
		os << (" Splitting results:       No") << endl;
		os << (" Result file:             ") << getOutput_file() << endl;
	}

}

void wig_cmd_option::verify() {
	std::ostringstream oss;
	oss << "Specified ";
	if (_format != "bed" && _format != "bowtie" && _format != "eland"
			&& _format != "sam" && _format != "bam") {
		oss << "file format " << _format
				<< " is not a valid entry. Please choose from : bowtie, sam, bam and bed."
				<< endl;
		throw not_in_range(oss.str().c_str());
	}

}

bool wig_cmd_option::isSplit() const {
	return _split;
}

void wig_cmd_option::setSplit(bool _split) {
	this->_split = _split;
}

bool wig_cmd_option::isGz() const {
	return _gz;
}

bool wig_cmd_option::isStranded() const {
	return _stranded;
}

void wig_cmd_option::setGz(bool _gz) {
	this->_gz = _gz;
}

void wig_cmd_option::setStranded(bool _stranded) {
	this->_stranded = _stranded;
}

uint32_t wig_cmd_option::getOverlapSz() const {
	return _overlap_sz;
}

uint32_t wig_cmd_option::getWinwdowSz() const {
	return _window_sz;
}

void wig_cmd_option::setOverlapSz(uint32_t overlapSz) {
	_overlap_sz = overlapSz;
}

void wig_cmd_option::setWinwdowSz(uint32_t winwdowSz) {
	_window_sz = winwdowSz;
}
