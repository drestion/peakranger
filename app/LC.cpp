/*
 * LC.cpp
 *
 *  Created on: Jul 10, 2012
 *      Author: xfeng
 */

#include "app/LC.h"
#include "common/stl_header.h"
#include "common/boost_header.h"
#include "common/ranger_debug.h"
#include "option_parser/LCOption.h"
#include "option_parser/OptionAux.h"
#include "utils/Tracer.h"
#include "bam_app/LibraryComplexity.h"
#include "bam_app/StockBamMultipleDatasetsApp.h"
#include "paired_reads_parser/BamFileSortOrderDetector.h"
#include "paired_reads_parser/BamBlockPEReadsParser2LImp.h"
#include "paired_reads_parser/BamBlockPEReadsParser1LImp.h"
#include "paired_reads_parser/PEBamFormatDetector.h"
//#include "utils/signaldumper.h"
using namespace std;
using namespace bam_app;
using namespace bam_app::aux;
using namespace boost;
using namespace options::aux;
using namespace parser;
using namespace parser::aux;
namespace app {
string LC::version = "";
LC::LC() {

}

LC::~LC() {

}

void LC::run(int argc, char** argv) {

	options::LCOption opts(version);
	try {
		opts.parse(argc, argv);
	} catch (boost::program_options::multiple_occurrences& e) {
		cout << "\nSome options were specified more than once.\n";
		exit(1);
	} catch (std::exception& e) {
		cout << "\n" << e.what() << "\n";
		cout << "\n" << "Provided args:\n" << printRawOpts(argc, argv) << "\n";
		exit(0);
	} catch (FileNotGood& e) {
		e.debugPrint();
		cout << "\n" << "Provided args:\n" << printRawOpts(argc, argv) << "\n";
		exit(0);
	} catch (...) {
		cout << "\n"
				<< "Unknown fatal exception caught while parsing command lines\n";
		cout << "\n" << "Provided args:\n" << printRawOpts(argc, argv) << "\n";
		exit(0);
	}

	utils::TimeStampTracer tracer(std::cout, opts.mVM.count("verbose"));
	string input(opts.mVM["data"].as<string>());
	vector<string> files(1, input);
	tracer << opts.printParsedOpts() << "\n";

	LibraryComplexity w(tracer);

	BamFileSortOrderDetector format;
	BamBlockPEReadsParserImp* parserImp;
	BamBlockPEReadsParser2LImp l2; // Could just use non-blocked parser.
	BamBlockPEReadsParser1LImp l1;
	PEBamFormatDetector pe;
	if (pe.isPairEndBamFile(input, 50000)) {
		uint32_t order = format.getType(input, 2);
		if (order & readname_sorted) {
			parserImp = &l2;
		} else {
			cout << "\nThe bam file is pair end and ";
			cout << "thus must be sorted by read names first.";
			cout << "\nTry:";
			cout << " \nsamtools sort -n " + input + " " + input;
			cout << ".sorted\n\n";
			exit(0);
		}
	} else {
		parserImp = &l1;
		w.mIsPE = false;
	}
	StockBamMultipleDatasetsApp lc(tracer);

	lc.setImp(&w);
	lc.setParserImp(parserImp);
	lc.run(files, cout);
	cout << "\n";

}

} /* namespace app */
