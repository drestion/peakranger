/*
 * ranger4ccmdparser.cpp
 *
 *  Created on: Feb 7, 2012
 *      Author: xfeng
 */

#include "ranger4ccmdparser.h"
#include "utils/stringutil.h"

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
#include "option_parser/OptionAux.h"
using namespace std;
using namespace boost;
using namespace utils;
using namespace options::aux;
using namespace boost::program_options;

namespace po = boost::program_options;

void ranger4c_cmd_parser::printHelp() const {
	cout << "Contact:\tXin Feng peak.ranger@gmail.com" << endl;
	//todo: ranger version
	cout << "\nChIP4C Version:" << version << endl;
	printf("\nSample usage:\n");
	printf("  ranger --format bowtie sample_file ");
	printf(" control_file -o result_file -t 4 --verbose\n");
	printf("\nData files:");
	printf("\n  -d,--data         data file.(REQUIRED) ");
//    printf("\n  -c,--control      control(input) file.(REQUIRED) ");
	printf("\n  --format          the format of the data file,");
	printf(
			"\n                    can be one of : bowtie, sam, bam and bed.(REQUIRED)");
	printf(
			"\n  --chr_table       process chromosomes contained in the specified chr table file. ");
	printf(
			"\n  --config          specify the location of the configuration file.");
	printf("\nQualities:");

	printf("\n  -p,--pval         p value cut off.(default:1e-4)");
	printf("\n  -q,--FDR          FDR cut off.(default:5e-2)");
	printf("\n  -l,--ext_length   read extension length.(default:200)");
	printf(
			"\n  -w,--window       window length used in window signal statistics.(default:5000)");
	printf(
			"\n  -r,--delta        sensitivity of summits detector, must be in the region(0, 1).(default:0.8)");
	printf("\n  -b,--bandwidth    bandwidth.(default:99)");
	printf(
			"\n  --pad             pad read coverage to avoid false positive summits.(default:false)");
	printf("\nRunning modes:");
	printf(
			"\n  --mode            specify the running mode, can be one of : region, resolution");
	printf("\n                    (default:region)");
	printf("\n  -t,--thread       number of threads.(default: 1)");

	printf("\nOutput:");
	printf(
			"\n  -o,--output       specify the location of output files.(REQUIRED)");
	printf(
			"\n  --nowig           do not generate wiggle(.wig) files.(default:off)");
	printf(
			"\n  --split,-s        generate one wig file per chromosome.(default:off)");
	printf("\nOther:");
	printf("\n  --verbose         print application progress.");
	printf("\n");

	exit(0);
}

void ranger4c_cmd_parser::parse() {

	opt all("Usage");
	opt other("Other");
	opt input("Input");
	opt output("Output");
	opt qualities("Qualities");
	opt running_modes("Running modes");

	other.add_options()("help", "show the usage")("verbose", "show progress")(
			"version", "output the version number");

	input.add_options()("data,d", po::value<string>(&_treat_dir),
			"chipseq data/treatment file")("control,c",
			po::value<string>(&_control_dir), "control file")("format",
			po::value<string>(&_format),
			"the format of the data file, can be one of : "
					"bowtie, eland, sam, bam and bed")("chr_table",
			po::value<string>(&_chr_table_file),
			"the file that contains the chromosomes to be processed")("config",
			po::value<string>(&_config_file), "configuration file");

	output.add_options()("output,o", po::value<string>(&_output_dir),
			"the output location")("nowig", "do not generate wiggle files")(
			"split,s", "generate one wig file per chromosome");

	qualities.add_options()

	("w,window", po::value<uint32_t>(&_w)->default_value(5000), "window length")(
			"pval,p", po::value<double>(&_p_cut_off)->default_value(1e-4),
			"p value cut-off")("FDR,q",
			po::value<double>(&_fdr_cut_off)->default_value(0.05),
			"FDR cut-off")("ext_length,l",
			po::value<uint32_t>(&_ext_length)->default_value(100),
			"read extension length")("delta,r",
			po::value<double>(&_delta)->default_value(8e-1),
			"sensitivity of the summit detector")("bandwidth,b",
			po::value<uint32_t>(&_bandwidth)->default_value(99),
			"smoothing bandwidth")("pad",
			"pad read coverage profile to avoid false positive summmits");

	running_modes.add_options()("mode",
			po::value<string>(&_mode)->default_value("region"),
			"specify the running mode, can be one of : region, resolution")(
			"thread,t", po::value<uint32_t>(&_no_of_thread)->default_value(1),
			"number of running threads")("binsize",
			po::value<uint32_t>(&_binlength)->default_value(10000),
			"bin length");

	if (_ac < 2) {
		printHelp();
		exit(0);
	}
	_config_file = "";
	p_opt popt;
	popt.add("data", 1).add("control", 1).add("output", 1);

	all.add(other).add(input).add(output).add(qualities).add(running_modes);

	po::store(
			po::command_line_parser(_ac, _av).options(all).positional(popt).run(),
			vm);
	po::notify(vm);

	if (_config_file != "") {
		ifstream ifs(_config_file.c_str());
		if (!ifs) {
			string str;
			str += "can not open config file: ";
			str += _config_file;
			throw RangerException(str);
		} else {
			po::store(po::parse_config_file(ifs, all), vm);
			po::notify(vm);
		}
	}

	setFormat(to_lower_copy(trim_copy(getFormat())));

	if (vm.count("help")) {
		setHelpRequested(true);
		printHelp();
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
	if (vm.count("control")) {
		this->setHas4CControl(true);
	} else {
		this->setHas4CControl(false);
		_control_dir = _treat_dir;
	}
	require("data", vm);
	require("output", vm);
	require("format", vm);

	file_r_good(_treat_dir.c_str());
	if (isHas4CControl()) {
		file_r_good(_control_dir.c_str());
	}

	string dir, file, file_ext;

	stringutil::get_dir_file(_treat_dir, dir, file, file_ext);

	setTreat_file(_treat_dir);
	setTreat_dir(dir);
	setTreatfilename(file);

	stringutil::get_dir_file(_control_dir, dir, file, file_ext);
	if (isHas4CControl()) {
		setControl_file(_control_dir);
		setControl_dir(dir);
		setControlfilename(file);
	} else {

		setControl_dir(dir);
		setControlfilename("simulated_" + file);
		setControl_file(dir + "simulated_" + file + file_ext);
	}

	file_w_good(_output_dir.c_str());
	stringutil::get_dir_file(_output_dir, dir, file, file_ext);
	//todo: linux only
	setTreat_wig_file(dir + "/" + _treatfilename + ".wig");
	setControl_wig_file(dir + "/" + _controlfilename + ".wig");
	setOutput_file(_output_dir);
	setOutput_dir(dir);

	if (vm.count("chr_table")) {
		file_r_good(_chr_table_file.c_str());
		setUsing_chr_table(true);
		setChr_table_file(_chr_table_file);
		parseChrTable(_chr_table_file, _chrs_to_parse);
	}

	if (vm.count("nowig")) {
		setNowig(true);
	} else {
		setNowig(false);
	}

	if (vm.count("pad")) {
		setPad(true);
	} else {
		setPad(false);
	}

	boost::trim(_mode);
	boost::to_lower(_mode);

	if (_mode == "resolution") {
		setDelta(0.2);
	}

	verify();

}

void ranger4c_cmd_parser::print_option(ostream & os) {
	os << ("ChIP4C version:           ") << version << endl;
	os << ("Data files:\n");
	os << (" File format:             ") << getFormat() << endl;
	os << (" Sample file:             ") << getTreat_file() << endl;
	if (isHas4CControl()) {
		os << (" Control file:            ") << getControl_file() << endl;
		if (getUsing_chr_table()) {
			os << (" Chr table file:          ") << getChr_table_file() << endl;
		//	string prefix(" Chrs to process:         ");
		//	os << (" Chrs to process:         ")
			//		<< utils::vector_to_string<string>(_chrs_to_parse, " ") << endl;

		}
	} else {
		os << (" Control file:          "
				"  Using simulated reads.(-c not specified)\n");
	}
	os << ("Qualities:\n");
	os << (" P value cut off:         ") << getCut_off() << endl;
	os << (" FDR cut off:             ") << getFdrCutOff() << endl;
	os << (" Window stat size:        ") << getW() << endl;
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
	os << (" Detection mode:          ") << getMode() << endl;
	os << (" Number of threads:       ") << getNo_of_thread() << endl;
	os << ("Output:\n");

	os << (" Regions:                 ") << getOutput_file() + "_region.bed"
			<< endl;
	os << (" Summits:                 ") << getOutput_file() + "_summit.bed"
			<< endl;
	os << (" Details of regions:      ") << getOutput_file() + "_details"
			<< endl;
	if (!getNowig()) {
		if (isSplit()) {
			os << (" Splitting results:       Yes") << endl;
			os << (" Data wiggle file:        ") << getTreat_wig_file() << endl;
			os
					<< ("                          and other splitted files in this directory")
					<< endl;
			os << (" Control wiggle file:     ") << getControl_wig_file()
					<< endl;
			os
					<< ("                          and other splitted files in this directory")
					<< endl;
		} else {

			os << (" Data wiggle file:        ") << getTreat_wig_file() << endl;
			os << (" Control wiggle file:     ") << getControl_wig_file()
					<< endl;
		}
	} else {
		os << (" Data wiggle file:        DISABLED(--nowig)\n");
		os << (" Control wiggle file:     DISABLED(--nowig)\n");
	}
}

void ranger4c_cmd_parser::verify() {
	std::ostringstream oss;
	oss << "Specified ";
	if (!is_in_range(_delta, 0.0, 1.0)) {
		oss << "delta value " << _delta << " is not in the range : (0, 1)"
				<< endl;
		throw not_in_range(oss.str().c_str());
	}

	if (!is_in_range(_p_cut_off, 0.0, 1.0)) {
		oss << "p-value cut off " << _p_cut_off
				<< " is not in the range : (0, 1)" << endl;
		throw not_in_range(oss.str().c_str());
	}

	if (_mode != "region" && _mode != "resolution") {
		oss << "running mode " << _mode
				<< " is not a valid entry. Please use either \"region\" or \"resolution\"."
				<< endl;
		throw not_in_range(oss.str().c_str());
	}

	if (_format != "bed" && _format != "bowtie" && _format != "eland"
			&& _format != "sam" && _format != "bam") {
		oss << "file format " << _format
				<< " is not a valid entry. Please choose from : bowtie, sam, bam and bed."
				<< endl;
		throw not_in_range(oss.str().c_str());
	}
	uint32_t maxThreads = boost::thread::hardware_concurrency();
	if (maxThreads < getNo_of_thread()) {
		oss << "number of threads " << getNo_of_thread()
				<< " is not valid. The current system allows up to "
				<< maxThreads << " threads.\n";
		throw not_in_range(oss.str().c_str());
	}
	//todo: fix the multi-threading issue in chip4c.cpp
	if (maxThreads < getNo_of_thread()) {
		oss << "number of threads " << getNo_of_thread()
				<< " is not valid. The current system allows up to " << 1
				<< " threads.\n";
		throw not_in_range(oss.str().c_str());
	}

}

bool ranger4c_cmd_parser::isSplit() const {
	return _split;
}

void ranger4c_cmd_parser::setSplit(bool _split) {
	this->_split = _split;
}

void ranger4c_cmd_parser::print_option_file(ostream & os) const {
	os << ("#Ranger version:           ") << version << endl;
	os << ("#Data files:\n");
	os << ("# File format:             ") << getFormat() << endl;
	os << ("# Sample file:             ") << getTreat_file() << endl;
	if (isHas4CControl()) {
		os << ("# Control file:            ") << getControl_file() << endl;
		if (getUsing_chr_table()) {
			os << ("# Chr table file:          ") << getChr_table_file()
					<< endl;
//			string prefix(" Chrs to process:         ");
	//		os << ("# Chrs to process:         ")
		//			<< utils::vector_to_string<string>(_chrs_to_parse," " ) << endl;

		}
	} else {
		os << ("# Control file:          "
				"  Using simulated reads.(-c not specified)\n");
	}
	os << ("#Qualities:\n");
	os << ("# P value cut off:         ") << getCut_off() << endl;
	os << ("# FDR cut off:             ") << getFdrCutOff() << endl;
	os << ("# Read extension length:   ") << getExt_length() << endl;
	os << ("# Smoothing bandwidth:     ") << getBandwidth() << endl;
	os << ("# Delta:                   ") << getDelta() << endl;
	os << ("# Pad region profile:      ");
	if (getPad()) {
		os << "#Enabled" << endl;
	} else {
		os << "#Disabled\n";
	}
	os << ("#Running modes:\n");
	os << ("# Detection mode:          ") << getMode() << endl;
	os << ("# Number of threads:       ") << getNo_of_thread() << endl;
	os << ("#Output:\n");

	os << ("# Regions:                 ") << getOutput_file() + "_region.bed"
			<< endl;
	os << ("# Summits:                 ") << getOutput_file() + "_summit.bed"
			<< endl;
	os << ("# Details of regions:      ") << getOutput_file() + "_details"
			<< endl;
	if (!getNowig()) {
		if (isSplit()) {
			os << ("# Splitting results:       Yes") << endl;
			os << ("# Data wiggle file:        ") << getTreat_wig_file()
					<< endl;
			os
					<< ("#                          and other splitted files in this directory")
					<< endl;
			os << ("# Control wiggle file:     ") << getControl_wig_file()
					<< endl;
			os
					<< ("#                          and other splitted files in this directory")
					<< endl;
		} else {

			os << ("# Data wiggle file:        ") << getTreat_wig_file()
					<< endl;
			os << ("# Control wiggle file:     ") << getControl_wig_file()
					<< endl;
		}
	} else {
		os << ("# Data wiggle file:        DISABLED(--nowig)\n");
		os << ("# Control wiggle file:     DISABLED(--nowig)\n");
	}
}

uint32_t ranger4c_cmd_parser::getW() const {
	return _w;
}

void ranger4c_cmd_parser::setW(uint32_t _w) {
	this->_w = _w;
}

bool ranger4c_cmd_parser::isHas4CControl() const {
	return _has4CControl;
}

void ranger4c_cmd_parser::setHas4CControl(bool _has4CControl) {
	this->_has4CControl = _has4CControl;
}
