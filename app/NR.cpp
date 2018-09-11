/*
 * NR.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */
#include "common/boost_test.h"
#include "common/stl_header.h"
#include "common/ranger_debug.h"
#include "common/boost_header.h"
#include "short_reads/reads.h"
#include "short_reads/readstools.h"
#include "region_detector/NoiseEstimator.h"
#include "utils/Tracer.h"
#include "app/NR.h"
#include "option_parser/NROption.h"
#include "option_parser/OptionAux.h"
#include "option_parser/cmd_option_parser.h"
#include "parser/bowtieParser.h"
#include "parser/samParser.h"
#include "parser/bamParser.h"
#include "parser/bedParser.h"
#include "parser/readsParser.h"

using namespace std;
using namespace boost;
using namespace reads;
using namespace options;
using namespace options::aux;
namespace app {
string NR::version = "";
NR::NR() {

}

NR::~NR() {

}

void NR::run(int argc, char** argv) {
	NROption opts(NR::version);
	try {
		opts.parse(argc, argv);
	} catch (std::exception& e) {
		cout << "\n" << e.what() << "\n";
		cout << "\nProvided args:\n" << printRawOpts(argc, argv) << "\n";
		exit(0);
	} catch (FileNotGood& e) {
		e.debugPrint();
		cout << "\nProvided args:\n" << printRawOpts(argc, argv) << "\n";
		exit(0);
	} catch (...) {
		cout
				<< "\nUnknown fatal exception caught while parsing command lines\n";
		cout << "\nProvided args:\n" << printRawOpts(argc, argv) << "\n";
		exit(0);
	}

	utils::Tracer tracer(cout, opts.mVM.count("verbose"));
	string input(opts.mVM["data"].as<string>());
	string control(opts.mVM["control"].as<string>());
	tracer << opts.printParsedOpts() << "\n";
	NoiseEstimator nre;
	nre.setFragSize(opts.mVM["ext_length"].as<uint32_t>());
	Reads treads, creads;
	boost::shared_ptr<readsParser> parser;
	if (opts.mVM["format"].as<string>() == cmd_option_parser::format_bowtie) {
		parser = boost::make_shared<bowtieParser>();
	} else if (opts.mVM["format"].as<string>()  == cmd_option_parser::format_sam) {
		parser = boost::make_shared<samParser>();
	} else if (opts.mVM["format"].as<string>()  == cmd_option_parser::format_bed) {
		parser = boost::make_shared<bedParser>();
	} else if (opts.mVM["format"].as<string>()  == cmd_option_parser::format_bam) {
		parser = boost::make_shared<bamParser>();
	} else {
		string str("The specified format ");
		str += opts.mVM["format"].as<string>() ;
		str += " has not been implemented yet.\n";
		throw std::logic_error(str.c_str());
	}

	parser->parse(treads, input);
	parser->parse(creads, control);
	reads_tools::verify_and_correct_Reads_both_strands(treads, creads);
	tracer << "Reads statistics:\n";
	tracer << " Treatment reads +:       " << treads.pos_reads.size() << "\n";
	tracer << " Treatment reads -:       " << treads.neg_reads.size() << "\n";
	tracer << " Average read length:     " << treads.getReadlength() << "\n";
	tracer << " Control reads +:         " << creads.pos_reads.size() << "\n";
	tracer << " Control reads -:         " << creads.neg_reads.size() << "\n";
	tracer << " Average read length:     " << creads.getReadlength() << "\n";
	double rate = nre.estimate(treads, creads);
	cout << "\nEstimated noise rate:" << rate << "\n\n";
}

} /* namespace app */
