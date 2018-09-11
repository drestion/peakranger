/*
 * chip4cwigcmdoption.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#include "chip4cwigcmdoption.h"
#include "utils/stringutil.h"
#include "option_parser/OptionAux.h"
#include "utils/exceptions.h"
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
using namespace std;
using namespace boost;
using namespace utils;
using namespace options::aux;
using namespace boost::program_options;


void chip4c_wig_cmd_option::printHelp() const {

    cout << "\nWig Version:" << version << endl;
    printf("\nSample usage:\n");
    printf("  wig --format=bowtie -d sample_file -o output_file");
    printf("\nData files:");
    printf("\n  -d                  data file.(REQUIRED) ");
    printf("\n  --format            the format of the data file,");
    printf("\n                      can be one of : bowtie, eland, sam, bam, bed and other.(REQUIRED)");
    printf("\n  --chr_table         only process chromosomes contained in the specified chr table file.");
    printf("\nQualities:");
    printf("\n  -l,--ext_length     read extension length.(default:200)");
    printf("\n  -w,                 window size.(default:10000)");
    printf("\n  -u,                 overlapping window size.(default:9000)");
    printf("\n  -r,--repeat         repeat zoomming times. (default:1)");
    printf("\nOutput:");
    printf("\n  -o                  specify the location of output files.(REQUIRED)");
    printf("\n  --split,-s          generate one wig file per chromosome.(default:off)");
    printf("\n  --strand,-x         generate one wig file per strand.(default:off)");
    printf("\nOther:");
    printf("\n  --verbose           print application progress.");
    printf("\n");

    exit(0);
}

void chip4c_wig_cmd_option::parse() {
    if (_ac < 2) {
        printHelp();
        exit(0);
    }


    opt all("Usage");
    opt other("Other");
    opt input("Input");
    opt output("Output");
    opt qualities("Qualities");

    other.add_options()

    ("help,h",
     "show the usage")

    ("verbose",
     "show progress")

    ("version",
     "output the version number");

    input.add_options()

    ("data,d",
     po::value<string>(&_treat_dir),
     "chipseq data/treatment file")

    ("format",
     po::value<string>(&_format),
     "the format of the data file, can be one of : "
     "bowtie, eland, sam, bam and bed")

    ("chr_table",
     po::value<string>(&_chr_table_file),
     "the file that contains the chromosomes to be processed");

    output.add_options()

    ("output,o",
     po::value<string>(&_output_dir),
     "the output location")

    ("split,s",
     "generate one wig file per chromosome")

    ("gzip,z",
     "gzip the result file")

    ("strand,x",
     "generate one wig file per strand");

    qualities.add_options()

    ("ext_length,l",
     po::value<uint32_t>(&_ext_length)->default_value(200),
     "read extension length")

    ("window_size,w",
     po::value<uint32_t>(&_window_sz)->default_value(10000),
     "window size zoom")

    ("window_overlap_size,u",
     po::value<uint32_t>(&_overlap_sz)->default_value((10000 - 1000)),
     "window size zoom overlap")

(     "repeat,r",
     po::value < uint32_t > (&_repeat)->default_value(5),
     "repeated zooming");

    p_opt popt;
    popt.add("data",
             1).add("output",
                    1);
    all.add(other).add(input).add(output).add(qualities);

    po::store(po::command_line_parser(_ac,
                                      _av).options(all).positional(popt).run(),
              vm);
    po::notify(vm);

    setFormat(to_lower_copy(trim_copy(getFormat())));

    if (vm.count("help")) {
        setHelpRequested(true);
        printHelp();
    }
    if (vm.count("verbose")) {
        setVerboseRequested(true);
    }

    if (vm.count("version")) {
        printVersion(cout,version);
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
    require("data",
            vm);
    require("output",
            vm);
    require("format",
            vm);

    file_r_good(_treat_dir.c_str());
    string dir, file, file_ext;

    stringutil::get_dir_file(_treat_dir,
                             dir,
                             file,
                             file_ext);

    setTreat_file(_treat_dir);
    setTreat_dir(dir);
    setTreatfilename(file);
    file_w_good(_output_dir.c_str());
    stringutil::get_dir_file(_output_dir,
                             dir,
                             file,
                             file_ext);
    //todo: linux only
    if (file_ext != "wig") {
        _output_dir += ".wig";
    }
    setOutput_file(_output_dir);
    setOutput_dir(dir);

    if (vm.count("chr_table")) {
        file_r_good(_chr_table_file.c_str());
        setUsing_chr_table(true);
        setChr_table_file(_chr_table_file);
        parseChrTable(_chr_table_file,
                      _chrs_to_parse);
    }
    verify();

}

void chip4c_wig_cmd_option::print_option_file(ostream & os) const {
    os << ("#Wig version:           ") << version << endl;
    os << ("Data files:\n");
    os << (" File format:             ") << getFormat() << endl;
    os << (" Sample file:             ") << getTreat_file() << endl;

    if (getUsing_chr_table()) {
        os << (" Chr table file:          ") << getChr_table_file() << endl;
      //  string prefix(" Chrs to process:         ");
      //  os << (" Chrs to process:         ")
      //  << utils::vector_to_string<string>(_chrs_to_parse, " ") << endl;

    }
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
void chip4c_wig_cmd_option::print_option(ostream & os) {

    os << ("Data files:\n");
    os << (" File format:             ") << getFormat() << endl;
    os << (" Sample file:             ") << getTreat_file() << endl;

    if (getUsing_chr_table()) {
        os << (" Chr table file:          ") << getChr_table_file() << endl;
        string prefix(" Chrs to process:         ");
        os << (" Chrs to process:         ")
        << vector_to_string<string>(_chrs_to_parse,
                                    (uint32_t) 3,
                                    prefix) << endl;

    }
    os << ("Qualities:\n");

    os << (" Read extension length:   ") << getExt_length() << endl;
    os << (" Window size:             ") << getWinwdowSz() << endl;
    os << (" Window overlap size:     ") << getOverlapSz() << endl;
    os << (" Repeat zooming:          ") << getRepeat()  << endl;

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

void chip4c_wig_cmd_option::verify() {

    std::ostringstream oss;
    oss << "Specified ";
    if (_format != "bed" && _format != "bowtie" && _format != "eland"
    && _format != "sam" && _format != "bam") {
        oss
        << "file format "
        << _format
        << " is not a valid entry. Please choose from : bowtie, eland, sam, bam and bed."
        << endl;
        throw not_in_range(oss.str().c_str());
    }
    oss.clear();
    oss.str("Specified ");
    if (_window_sz <= _overlap_sz) {
        oss <<
        " window size " << _window_sz << " is smaller than the window"
        << " overlap size " << _overlap_sz << "\n";
        throw not_in_range(oss.str().c_str());
    }

}

bool chip4c_wig_cmd_option::isSplit() const {
    return _split;
}

void chip4c_wig_cmd_option::setSplit(bool _split) {
    this->_split = _split;
}

bool chip4c_wig_cmd_option::isGz() const {
    return _gz;
}

bool chip4c_wig_cmd_option::isStranded() const {
    return _stranded;
}

void chip4c_wig_cmd_option::setGz(bool _gz) {
    this->_gz = _gz;
}

void chip4c_wig_cmd_option::setStranded(bool _stranded) {
    this->_stranded = _stranded;
}

uint32_t chip4c_wig_cmd_option::getOverlapSz() const {
    return _overlap_sz;
}

uint32_t chip4c_wig_cmd_option::getWinwdowSz() const {
    return _window_sz;
}

void chip4c_wig_cmd_option::setOverlapSz(uint32_t overlapSz) {
    _overlap_sz = overlapSz;
}

void chip4c_wig_cmd_option::setWinwdowSz(uint32_t winwdowSz) {
    _window_sz = winwdowSz;
}

uint32_t chip4c_wig_cmd_option::getRepeat() const
{
    return _repeat;
}

void chip4c_wig_cmd_option::setRepeat(uint32_t _repeat)
                                      {
    this->_repeat = _repeat;
}

inline void chip4c_wig_cmd_option::require(const char* opt,
                                           po::variables_map& vm) {
    if (!vm.count(opt)) {
        string _opt(opt);
        string oopt;
        if (_opt.size() > 1) {
            oopt += "--";
        } else {
            oopt += "-";
        }
        oopt += _opt;
        cerr << "Missing required option : " << oopt << endl;
        exit(1);
    }

}
