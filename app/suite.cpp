/*
 * suite.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "utils/exceptions.h"
#include "app/PeakRanger.h"
#include "app/CCAT.h"
#include "app/NR.h"
#include "app/WigPE.h"
#include "app/Wig.h"
#include "app/LC.h"
#include "app/bcp.h"

using namespace std;
using namespace app;

void helpAndQuit(ostream& os) {
	os << "\n";
	os << "Program: peakranger (Resources for NGS data)\n";
	os << "Version: " << VERSION << "\n";
	os << "Contact: Xin Feng (peak.ranger@gmail.com)\n\n";
	os << "Usage:   peakranger <command> <arguments>\n\n";
	os << "Commands:                                   \n";
	os << "         nr      estimate data quality\n";
	os << "         lc      calculate library complexity\n";
	os << "         wig     generate wiggle files\n";
	os << "         wigpe   generate wiggle files for paired reads\n";
	os << "         ranger  peak calling for sharp peaks\n";
	os << "         ccat    peak calling for broad peaks\n";
	os << "         bcp     peak calling for complex broad peaks\n";
	os << "\n";
	exit(0);
}

int main(int argc, char **argv) {
	if (argc < 2){
		helpAndQuit(cout);
	}
	stringstream version;
	version << VERSION;

	string cmd(argv[1]);
	boost::trim(cmd);
	boost::to_lower(cmd);

	NR::version = version.str();
	Wig::version = version.str();
	WigPE::version = version.str();
	PeakRanger::version = version.str();
	CCAT::version = version.str();
	LC::version = version.str();
	BCP::version = version.str();

	try {
		if (cmd == "nr")
			NR::run(argc - 1, argv + 1);
		else if (cmd == "wig")
			Wig::run(argc - 1, argv + 1);
		else if (cmd == "wigpe")
			WigPE::run(argc - 1, argv + 1);
		else if (cmd == "ranger")
			PeakRanger::run(argc - 1, argv + 1);
		else if (cmd == "ccat")
			CCAT::run(argc - 1, argv + 1);
		else if (cmd == "lc")
			LC::run(argc - 1, argv + 1);
		else if (cmd == "bcp")
			BCP::run(argc - 1, argv + 1);
		else {
			cout << "\nUnrecognized command name: " << cmd;
			helpAndQuit(cout);
		}
	} catch (boost::program_options::multiple_occurrences& e) {
		cout << "\nSome options were specified more than once.\n";
		exit(1);
	} catch (std::bad_alloc &e) {
		cout << "\nNot enough memory. std::bad_alloc caught\n";
		exit(1);
	} catch (RangerException &e) {
		cout << "\n";
		e.debugPrint();
		cout << "\n";
		exit(1);
	} catch (std::exception &e) {
		cout << "\n" << e.what() << "\n";
		exit(1);
	} catch (...) {
		cout << "\nUnknown error" << endl;
		exit(1);
	}
	return 0;
}
