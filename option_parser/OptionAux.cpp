/*
 * OptionAux.cpp
 *
 *  Created on: Jun 25, 2012
 *      Author: xfeng
 */

#include "option_parser/OptionAux.h"
#include "common/ranger_debug.h"
#include "utils/Guarded.h"
using namespace std;
using namespace boost::program_options;
namespace options {
namespace aux {

opt SharedOpts::help_verbose_version("Help");

SharedOpts::SharedOpts() {
	help_verbose_version.add_options()("help,h", "show help info")("verbose",
			"show progress when possible")("version",
			"output the version number");
}

void printHelp(const options_description& opts, std::ostream& os) {
	os << opts;
}

void file_r_good(const char* file) {
	ifstream ofs(file);
	utils::Guarded<FileNotGood> fpg(!ofs, file);
}

void require(const char* opt, variables_map& vm) {
	if (!vm.count(opt)) {
		string _opt(opt);
		string oopt;
		if (_opt.size() > 1) {
			oopt += "--";
		} else {
			oopt += "-";
		}
		oopt += _opt;
		throw std::logic_error("Missing required option : " + oopt);
	}
}

void file_w_good(const char* file) {
	ofstream ofs(file, ios_base::out);
	utils::Guarded<FileNotGood> fpg(!ofs, file);
	remove(file);
}

void printVersion(ostream& os, const string& version) {
	os << "\n";
	os << version;
	os << "\n\n";
	exit(0);
}
void parseChrTable(string& filename, vector<string>& results) {
	ifstream ifs(filename.c_str());

	if (!(ifs.good())) {
		string file("Chromosome table file ");
		file += filename;
		throw FileNotGood(file.c_str());
	}
	string line;
	results.resize(0);
	while (getline(ifs, line)) {
		boost::trim(line);
		boost::to_lower(line);
		results.push_back(line);
	}
}

void requireNotEqual(const char* opt, const char* notequal,
		const char* exepctedvalue, boost::program_options::variables_map& vm) {
	if (vm.count(opt)) {
		if (vm[opt].as<string>() == string(notequal)) {
			string actual(vm[opt].as<string>());
			string _opt(opt);
			string oopt;
			if (_opt.size() > 1) {
				oopt += "--";
			} else {
				oopt += "-";
			}
			oopt += _opt;
			throw std::logic_error(
					"The option : " + oopt + " has an illegal value:" + actual
							+ ". Expected value:" + string(exepctedvalue));
		}
	}
}
std::string printRawOpts(int argc, char** argv) {
	stringstream ss;
	for (int i = 0; i < argc; i++) {
		ss << argv[i];
		if (i + 1 < argc) {
			ss << " ";
		}
	}
	return ss.str();
}

}
}

/* namespace options */
