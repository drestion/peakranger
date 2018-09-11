/*
 * Wig.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */
#include <stdexcept>
#include <fstream>
#include "parser/readsParser.h"
#include "utils/exceptions.h"
#include "parser/bowtieParser.h"
#include "parser/samParser.h"
#include "parser/bamParser.h"
#include "parser/bedParser.h"
#include "short_reads/reads.h"
#include "option_parser/wig_cmd_option.h"
#include "result_reporter/result_reporter.h"
#include "result_reporter/bed6_result_reporter.h"
#include "wiggle/wiggle_reporter.h"
#include "wiggle/JTwigglefile.h"
#include "wiggle/strandedjtwiggle.h"
#include "utils/logger.h"
#include "utils/assert_helpers.h"
#include "utils/Tracer.h"
#include "common/boost_header.h"
#include "app/Wig.h"
using namespace std;
using namespace boost;

namespace app {
std::string Wig::version = "";
Wig::Wig() {

}

Wig::~Wig() {

}

void Wig::run(int argc, char** argv) {
#ifndef DEBUG
	try {
#endif

		wig_cmd_option option(argc, argv, version);
		option.parse();
		utils::Tracer tracer(cout, option.getVerboseRequested());
		if (option.getVerboseRequested()) {
			option.print_option(cout);
		}
		boost::shared_ptr<readsParser> parser;

		if (option.getFormat() == cmd_option_parser::format_bowtie) {
			parser = boost::make_shared<bowtieParser>();
		} else if (option.getFormat() == cmd_option_parser::format_sam) {
			parser = boost::make_shared<samParser>();
		} else if (option.getFormat() == cmd_option_parser::format_bed) {
			parser = boost::make_shared<bedParser>();
		} else if (option.getFormat() == cmd_option_parser::format_bam) {
			parser = boost::make_shared<bamParser>();
		} else {
			string str("The specified format ");
			str += option.getFormat();
			str += " has not been implemented yet.\n";
			throw not_in_range(str.c_str());
		}

		Reads reads;
		string ga = option.getTreat_file();
		if (option.getChrtableSpecified()) {
			vector<string> chrs_to_parse = option.getChrs_to_parse();
			parser->parse(reads, ga, chrs_to_parse);
		} else {
			parser->parse(reads, ga);
		}
		tracer << "\nReads statistics:\n";
		tracer << "\n Processed reads +:       " << reads.pos_reads.size();
		tracer << "\n Processed reads -:       " << reads.neg_reads.size();
		tracer << "\n Average read length:     " << reads.getReadlength()
				<< "\n";
		boost::shared_ptr<wiggle_reporter> wig;
		if (option.isStranded()) {
			wig = boost::make_shared<stranded_jtwiggle>();
		} else {
			wig = boost::make_shared<JT_wiggle_file>();
		}
		wig->use_default_setting();
		wig->setReadextlength(option.getExt_length());
		if (option.getExt_length() < reads.getReadlength()) {
			tracer << "Warning: Specified read extension length "
					<< option.getExt_length() << " is shorter than"
					<< " the read length " << reads.getReadlength()
					<< ". Forced to use the read length"
					<< " as the read extension length\n ";
			wig->setReadextlength(reads.getReadlength());
		}
		wig->setReadlength(reads.getReadlength());
		ga = option.getOutput_file();
		if (option.isSplit()) {
			if (option.isGz()) {
				wig->split_export_wiggle_gzip(reads, ga.c_str());
			} else {
				wig->split_export_wiggle(reads, ga.c_str());
			}
		} else {
			if (option.isGz()) {
				wig->export_wiggle_gzip(reads, ga.c_str());
			} else {
				wig->export_wiggle(reads, ga.c_str());
			}
		}
		tracer << "\n\nWiggle file generation completed.\n\n";

#ifndef DEBUG
	} catch (boost::program_options::multiple_occurrences& e) {
		cout << "\nSome options were specified more than once.\n";
		exit(1);
	} catch (std::exception &e) {
		cout << "\n" << e.what() << "\n";
		exit(1);
	} catch (RangerException &e) {
		cout << "\n";
		e.debugPrint();
		cout << "\n";
		exit(1);
	} catch (...) {
		cout << "Unknown error" << endl;
		exit(1);
	}
#endif
}

} /* namespace app */
