/*
 * CCAT.cpp
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
#include "option_parser/cmd_option_parser.h"
#include "region_detector/region_detector.h"
#include "region_detector/calledpeak.h"
#include "region_detector/ccat.h"
#include "result_reporter/result_reporter.h"
#include "result_reporter/bed6_result_reporter.h"
#include "wiggle/wiggle_reporter.h"
#include "wiggle/JTwigglefile.h"
#include "utils/logger.h"
#include "utils/stringutil.h"
#include "utils/Stamp.h"
#include "utils/timer.h"
#include "short_reads/readstools.h"
#include "ggplay/chipseqhtmlreporter.h"
#include "app/CCAT.h"
#include "option_parser/CCATOption.h"
#include "option_parser/OptionAux.h"
#include "option_parser/ccatcmdoptionparser.h"
#include "app/AppAux.h"
#include "tab_file/NearbyGeneFinder.h"
#include "concepts/RegionInt32.h"
#include "common/boost_header.h"
using namespace tab_file;
using namespace std;
using namespace boost;
using namespace utils;
using namespace options;
using namespace options::aux;
using namespace ranger::concepts;
#define foreach BOOST_FOREACH

typedef map<string, vector<called_peak> > enriched_regions;
typedef vector<called_peak>::iterator ritrr;
typedef enriched_regions::iterator pritrr;

namespace {
NearbyGeneFinder _nbgf;
}
namespace app {
std::string CCAT::version = "";
void CCAT::run(int argc, char** argv) {
#ifndef DEBUG
	try {
#endif

		options::ccat_cmd_option_parser option(argc, argv, version);

		try {
			option.parse();
		} catch (std::exception& e) {
			cout << "\n" << e.what() << "\n";
			cout << "\n" << "Provided args:\n" << printRawOpts(argc, argv)
					<< "\n";
			exit(0);
		} catch (FileNotGood& e) {
			e.debugPrint();
			cout << "\n" << "Provided args:\n" << printRawOpts(argc, argv)
					<< "\n";
			exit(0);
		}

		SET_LOG_FILE(option.getOutput_file()+"_ccat.log");

		SET_LOG_LEVEL("DEBUG3");

		LOG_DEBUG1("\n\n************Starting an new round***********\n\n");

		_nbgf.setAnnoFile(option.getGeneAnnoFile());
		_nbgf.setSearchSpan(option.getHtmlRegionLength());
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

		boost::shared_ptr<region_detector> detector =
				boost::make_shared<ccat>();
		boost::shared_ptr<result_reporter> reporter = boost::make_shared<
				bed6_result_reporter>();
		Reads treads, creads;
		string ga = option.getTreat_file();
		parser->parse(treads, ga);
		tracer << "\nReads statistics:\n";
		tracer << "\n Treatment reads +:       " << treads.pos_reads.size();
		tracer << "\n Treatment reads -:       " << treads.neg_reads.size();
		tracer << "\n Average read length:     " << treads.getReadlength();
		ga = option.getControl_file();
		parser->parse(creads, ga);
		tracer << "\n Control reads +:         " << creads.pos_reads.size();
		tracer << "\n Control reads -:         " << creads.neg_reads.size();
		tracer << "\n Average read length:     " << creads.getReadlength();
		tracer << "\n Verifying reads...\n";
		reads_tools::verify_and_correct_Reads_both_strands(treads, creads);
		tracer << "\nReads statistics after correction:\n";
		tracer << "\n Treatment reads +:       " << treads.pos_reads.size();
		tracer << "\n Treatment reads -:       " << treads.neg_reads.size();
		tracer << "\n Control reads +:         " << creads.pos_reads.size();
		tracer << "\n Control reads -:         " << creads.neg_reads.size();
		if (treads.size() < 1) {
			tracer << "\n\nNo reads were found in the treatment data. "
					<< "\nccat must stop here.\n\n";
			exit(0);
		}
		if (creads.size() < 1) {
			tracer << "\n\nNo reads were found in the control data. "
					<< "\nccat must stop here.\n\n";
			exit(0);
		}
		uint32_t treadslen = treads.getReadlength();
		uint32_t creadslen = creads.getReadlength();
		uint32_t optExtlen = option.getExt_length();
		uint32_t longerlen = treadslen > creadslen ? treadslen : creadslen;
		if (optExtlen < treadslen || optExtlen < creadslen) {
			tracer << "\nWarning: Specified read extension length " << optExtlen
					<< " is shorter than" << " the read length " << treadslen
					<< "(+)/" << creadslen << "(-). Forced to use " << longerlen
					<< " as the read extension length\n ";
			option.setExt_length(longerlen);
		}
		{
			tracer << "\n\n Calling peaks...\n\n";
			ga = option.getOutput_file() + "_raw";
			detector->detectSummits(treads, creads, option);
		}

		size_t fdr_passed = 0;
		size_t fdr_failed = 0;
		{
			ga = (option.getOutput_file() + "_region.bed");
			ofstream of(ga.c_str());
			if (!(of.is_open())) {
				throw FileNotGood(ga.c_str());
			}
			ga = (option.getOutput_file() + "_summit.bed");
			ofstream of_smt(ga.c_str());
			if (!(of_smt.is_open())) {
				throw FileNotGood(ga.c_str());
			}
			ga = (option.getOutput_file() + "_details");
			ofstream of_raw(ga.c_str());
			if (!(of_raw.is_open())) {
				throw FileNotGood(ga.c_str());
			}
			utils::Stamp::citationRangerCCATAndDate(of_raw);
			utils::Stamp::citationRangerCCATAndDate(of_smt);
			utils::Stamp::citationRangerCCATAndDate(of);
			option.print_option_file(of_raw);
			option.print_option_file(of_smt);
			option.print_option_file(of);

			of_smt << "\n\n#summit_chr\tsummit_start\tsummit_end"
					"\tsummit_ID\tsummit_FDR\tsummit_strand\n";

			of << "\n\n#region_chr\tregion_start\tregion_end"
					"\tregion_ID\tregion_fdr\tregion_strand\n";

			of_raw << "\n#region_chr\tregion_start\t" << "region_end\t"
					<< "nearby_genes(" << option.getHtmlRegionLength() / 1000
					<< "kbp)" << "\tregion_ID\t" << "region_summits\t"
					<< "region_fdr\tregion_strand\tregion_treads\tregion_creads\n";

			pritrr it = detector->_resultRegions.begin();
			for (; it != detector->_resultRegions.end(); it++) {

				foreach(called_peak pk , it->second) {
					of_raw << it->first << "\t" << pk.first << "\t" << pk.second
							<< "\t";
					if (option.needHtml()) {
						std::stringstream geneNamess;
						vector<TabGene> genes;
						_nbgf.getOverlappedGenes(it->first,
								RegionUint32(pk.first, pk.second), genes);
						foreach(TabGene& g, genes) {
							geneNamess << g.name << ",";
						}
						string geneNames(geneNamess.str());
						boost::replace_last(geneNames, ",", "");
						of_raw << geneNames;
					}
					of_raw << "\tccat";
					if (pk.q <= option.getFdrCutOff()) {
						of_raw << "_fdrPassed_" << fdr_passed;
					} else {
						of_raw << "_fdrFailed_" << fdr_failed;
					}
					of_raw << "_fdr_" << pk.q;
					of_raw << "\t" << utils::vector_to_string(pk.summits, ",");
					vector<uint32_t>::iterator sit = pk.summits.begin();
					for (; sit != pk.summits.end(); sit++) {

						of_smt << it->first << "\t" << *sit << "\t"
								<< (*sit) + 1 << "\tccat_region_" << pk.first
								<< "_" << pk.second;
						if (pk.q <= option.getFdrCutOff()) {
							of_smt << "_fdrPassed_" << fdr_passed;
						} else {
							of_smt << "_fdrFailed_" << fdr_failed;
						}

						of_smt << "\t" << pk.q << "\t+\n";
					}

					of_raw << "\t" << pk.q << "\t+\t" << pk.treads << "\t"
							<< pk.creads << "\n";

					of << it->first << "\t" << pk.first << "\t" << pk.second
							<< "\tccat";
					if (pk.q <= option.getFdrCutOff()) {
						of << "_fdrPassed_" << fdr_passed++;
					} else {
						of << "_fdrFailed_" << fdr_failed++;
					}
					of << "_fdr_" << pk.q << "\t" << pk.q << "\t+\n";
				}
			}
		}

		tracer << "\n\nTotal regions discovered:\t" << fdr_failed + fdr_passed;
		tracer << "\n\nTotal regions passed FDR cutoff:\t" << fdr_passed;

		if (option.needHtml()) {
			tracer << "\n\n Generating reports:\n";
			chipseq_html_reporter rptr;
			utilprint::ccatCitation ct;
			utilprint::citation ct2;
			rptr.addCitation(ct2.tostring() + " and " + ct.tostring());
			rptr.setRegionLength(option.getHtmlRegionLength());
			map<string, vector<called_peak> > passFDR;
			app::aux::filterByFDR(detector->_resultRegions, passFDR,
					option.getFdrCutOff());

			rptr.generate_report(treads, creads, passFDR, option);

		}
		tracer << "\n\nProgram finished.\n\n";

		LOG_DONE();

#ifndef DEBUG
	} catch (boost::program_options::multiple_occurrences& e) {
		cout << "\nSome options were specified more than once.\n";
		exit(1);
	} catch (std::bad_alloc &e) {
		cout << "\nNot enough memory. std::bad_alloc\n";
		exit(1);
	} catch (std::exception &e) {
		cout << "\n" << e.what();
		exit(1);
	} catch (RangerException &e) {
		cout << "\n";
		e.debugPrint();
		cout << "\n";
		exit(1);
	} catch (...) {
		cout << "\nUnknown error" << endl;
		exit(1);
	}
#endif

}

CCAT::~CCAT() {

}

} /* namespace app */

