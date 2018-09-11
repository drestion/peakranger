/*
 * PeakRanger.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#include "app/PeakRanger.h"
#include "app/AppAux.h"
#include <stdexcept>
#include <fstream>

#include "common/boost_header.h"
#include "parser/readsParser.h"
#include "utils/exceptions.h"
#include "parser/bowtieParser.h"
#include "parser/samParser.h"
#include "parser/bamParser.h"
#include "parser/bedParser.h"
#include "option_parser/cmd_option_parser.h"
#include "option_parser/peakranger_cmd_option_parser.h"
#include "region_detector/region_detector.h"
#include "region_detector/calledpeak.h"
#include "region_detector/fdr_based_thresholder.h"
#include "result_reporter/result_reporter.h"
#include "result_reporter/bed6_result_reporter.h"
#include "wiggle/wiggle_reporter.h"
#include "wiggle/JTwigglefile.h"
#include "tab_file/NearbyGeneFinder.h"
#include "utils/logger.h"
#include "utils/Stamp.h"
#include "utils/stringutil.h"
#include "utils/util_print.h"
#include "utils/Tracer.h"
#include "short_reads/readstools.h"
#include "concepts/RegionInt32.h"
#include "ggplay/chipseqhtmlreporter.h"
#include "region_detector/ccat_main.h"
using namespace std;
using namespace boost;
using namespace utils;
using namespace tab_file;
using namespace ranger::concepts;

#define foreach BOOST_FOREACH

typedef map<string, vector<called_peak> > enriched_regions;
typedef vector<called_peak>::iterator ritrr;
typedef enriched_regions::iterator pritrr;

namespace {
bool sorter_by_pval(called_peak r1, called_peak r2) {
	return r1.p < r2.p;
}

NearbyGeneFinder _nbgf;
}
namespace app {

PeakRanger::PeakRanger() {

}

PeakRanger::~PeakRanger() {

}

void CalculateFDR(peakranger_cmd_option_parser & option,
		boost::shared_ptr<region_detector> & detector) {
	if (option.getVerboseRequested()) {
		cout << "\n Calculating FDR...";
	}
	string chr;
	pritrr it = detector->_resultRegions.begin();
	for (; it != detector->_resultRegions.end(); it++) {
		chr = it->first;
		size_t d_pk_cnt = it->second.size();
		if (d_pk_cnt < 1) {
			continue;
		}
		std::sort(detector->_resultRegions[chr].begin(),
				detector->_resultRegions[chr].end(), sorter_by_pval);
		size_t j = 0;
		size_t _rk = 0;
		vector<called_peak> & _pk = detector->_resultRegions[chr];
		for (; j < _pk.size(); j++) {
			if (_pk[j].p < 0L) {
				_pk[j].q = 0;
			} else {
				_rk++;
				_pk[j].q = _pk[j].p * _pk.size() / _rk;
			}
		}

	}
	if (option.getVerboseRequested()) {
		cout << "finished\n";
	}

}

void OutputResults(peakranger_cmd_option_parser & option,
		boost::shared_ptr<region_detector> & detector, size_t & fdr_passed,
		size_t & fdr_failed) {
	string ga;
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
	utils::Stamp::citationAndDate(of_raw);
	utils::Stamp::citationAndDate(of_smt);
	utils::Stamp::citationAndDate(of);
	option.print_option_file(of_raw);
	option.print_option_file(of_smt);
	option.print_option_file(of);
	of_smt
			<< "\n#summit_chr\tsummit_start\tsummit_end\tsummit_ID\tsummit_FDR\tsummit_strand\n";
	of
			<< "\n#region_chr\tregion_start\tregion_end\tregion_ID\tregion_FDR\tregion_strand\n";
	of_raw << "\n#region_chr\tregion_start\t" << "region_end\t"
			<< "nearby_genes(" << option.getHtmlRegionLength() / 1000 << "kbp)"
			<< "\tregion_ID\t" << "region_summits\tregion_pvalue\t"
			<< "region_FDR\tregion_strand\tregion_treads\tregion_creads\n";

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
			of_raw << "\tranger";
			if (pk.q <= option.getFdrCutOff()) {
				of_raw << "_fdrPassed_" << fdr_passed;
			} else {
				of_raw << "_fdrFailed_" << fdr_failed;
			}
			of_raw << "_pval_" << pk.p << "_fdr_" << pk.q;
			of_raw << "\t" << utils::vector_to_string(pk.summits, ",");
			vector<uint32_t>::iterator sit = pk.summits.begin();
			for (; sit != pk.summits.end(); sit++) {

				of_smt << it->first << "\t" << *sit << "\t" << (*sit) + 1
						<< "\tranger_region_" << pk.first << "_" << pk.second
						<< "_pval_" << pk.p;
				if (pk.q <= option.getFdrCutOff()) {
					of_smt << "_fdrPassed_" << fdr_passed;
				} else {
					of_smt << "_fdrFailed_" << fdr_failed;
				}

				of_smt << "\t" << pk.q << "\t+\n";
			}

			of_raw << "\t" << pk.p << "\t" << pk.q << "\t+\t" << pk.treads
					<< "\t" << pk.creads << "\n";

			of << it->first << "\t" << pk.first << "\t" << pk.second
					<< "\tranger";
			if (pk.q <= option.getFdrCutOff()) {
				of << "_fdrPassed_" << fdr_passed++;
			} else {
				of << "_fdrFailed_" << fdr_failed++;
			}
			of << "_pval_" << pk.p << "_fdr_" << pk.q << "\t" << pk.q
					<< "\t+\n";
		}
	}
}

void getParser(peakranger_cmd_option_parser & option,
		boost::shared_ptr<readsParser> & parser) {
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

}

void parseReads(peakranger_cmd_option_parser& option, string readsfile,
		boost::shared_ptr<readsParser>& parser, Reads& treads) {

	if (option.getChrtableSpecified()) {
		vector<string> chrs_to_parse = option.getChrs_to_parse();
		parser->parse(treads, readsfile, chrs_to_parse);
	} else {
		parser->parse(treads, readsfile);
	}
}
std::string PeakRanger::version = "";
void PeakRanger::run(int argc, char** argv) {

#ifndef DEBUG
	try {
#endif

		peakranger_cmd_option_parser option(argc, argv, version);
		option.parse();
		SET_LOG_FILE(option.getOutput_file()+"_ranger.log");

		SET_LOG_LEVEL("DEBUG3");

		LOG_DEBUG1("\n\n************Starting an new round***********\n\n");

		_nbgf.setAnnoFile(option.getGeneAnnoFile());
		_nbgf.setSearchSpan(option.getHtmlRegionLength());
		utils::Tracer tracer(cout, option.getVerboseRequested());

		if (option.getVerboseRequested()) {
			option.print_option(cout);
		}

		boost::shared_ptr<readsParser> parser;
		getParser(option, parser);
		boost::shared_ptr<region_detector> detector = boost::make_shared<
				fdr_based_thresholder>();

		boost::shared_ptr<wiggle_reporter> wig = boost::make_shared<
				JT_wiggle_file>();
		wig->use_default_setting();

		Reads treads, creads;
		parseReads(option, option.getTreat_file(), parser, treads);

		tracer << "Reads statistics:\n";
		tracer << " Treatment reads +:       " << treads.pos_reads.size()
				<< "\n";
		tracer << " Treatment reads -:       " << treads.neg_reads.size()
				<< "\n";
		tracer << " Average read length:     " << treads.getReadlength()
				<< "\n";
		parseReads(option, option.getControl_file(), parser, creads);

		tracer << " Control reads +:         " << creads.pos_reads.size()
				<< "\n";
		tracer << " Control reads -:         " << creads.neg_reads.size()
				<< "\n";
		tracer << " Average read length:     " << creads.getReadlength()
				<< "\n";

		tracer << " Verifying reads...\n";
		reads_tools::verify_and_correct_Reads_both_strands(treads, creads);
		tracer << "Reads statistics after correction:\n";
		tracer << " Treatment reads +:       " << treads.pos_reads.size()
				<< "\n";
		tracer << " Treatment reads -:       " << treads.neg_reads.size()
				<< "\n";
		tracer << " Control reads +:         " << creads.pos_reads.size()
				<< "\n";
		tracer << " Control reads -:         " << creads.neg_reads.size()
				<< "\n";
		if (treads.size() < 1) {
			tracer << "\nNo reads were found in the treatment data. "
					"Ranger must stop here.\n\n";
			exit(0);
		}
		if (creads.size() < 1) {
			tracer << "\nNo reads were found in the control data. "
					"Ranger must stop here.\n\n";
			exit(0);
		}
		uint32_t treadslen = treads.getReadlength();
		uint32_t creadslen = creads.getReadlength();
		uint32_t optExtlen = option.getExt_length();
		uint32_t longerlen = treadslen > creadslen ? treadslen : creadslen;
		if (optExtlen < treadslen || optExtlen < creadslen) {
			tracer << "Warning: Specified read extension length " << optExtlen
					<< " is shorter than" << " the read length " << treadslen
					<< "(+)/" << creadslen << "(-). Forced to use " << longerlen
					<< " as the read extension length\n ";
			option.setExt_length(longerlen);
		}
		{
			tracer << "\n Calling peaks...\n" << "\n";
			string ga = (option.getOutput_file() + "_raw");
			detector->detectSummits(treads, creads, option);
		}
		CalculateFDR(option, detector);

		size_t fdr_passed = 0;
		size_t fdr_failed = 0;
		OutputResults(option, detector, fdr_passed, fdr_failed);
		tracer << "\nTotal regions discovered:\t" << fdr_failed + fdr_passed;
		tracer << "\nTotal regions passed FDR cutoff:\t" << fdr_passed << "\n";
		if (option.needHtml()) {
			tracer << "\n Generating reports:\n";
			chipseq_html_reporter rptr;
			utilprint::citation ct;
			rptr.addCitation(ct.tostring());
			rptr.setRegionLength(option.getHtmlRegionLength());
			map<string, vector<called_peak> > passFDR;
			app::aux::filterByFDR(detector->_resultRegions, passFDR,
					option.getFdrCutOff());
			rptr.generate_report(treads, creads, detector->_resultRegions,
					option);
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
		cout << "\n" << e.what() << "\n";
		exit(1);
	} catch (RangerException &e) {
		cout << "\n";
		e.debugPrint();
		cout << "\n";
		exit(1);
	} catch (...) {
		cout << "Unknown error" << "\n";
		exit(1);
	}
#endif
}

} /* namespace app */
