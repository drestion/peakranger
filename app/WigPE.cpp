/*
 * WigPE.cpp
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */
#include "common/stl_header.h"
#include "common/boost_header.h"
#include "common/ranger_debug.h"
#include "bam_app/OnlineBamMultiReportApp.h"
#include "bam_app/OnlineBamMultiReportAppImp.h"
#include "bam_app/WigPEOnline.h"
#include "option_parser/WigPEOption.h"
#include "option_parser/OptionAux.h"
#include "utils/Tracer.h"
#include "app/WigPE.h"
#include <boost/algorithm/string.hpp>
using namespace std;
using namespace bam_app;
using namespace boost;
using namespace options::aux;
namespace app {

WigPE::WigPE() {

}

WigPE::~WigPE() {

}
std::string WigPE::version = "";
void WigPE::run(int argc, char** argv) {

	options::WigPEOption opts(version);
	try {
		opts.parse(argc, argv);
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
	string output(normOutFileName(opts.mVM["output"].as<string>()));
	tracer << opts.printParsedOpts() << "\n";

	aux::WigPEOnline w(tracer);

	w.setExt(opts.mVM["ext_length"].as<uint32_t>());
	w.setSplitByChr(opts.mVM.count("split"));
	w.setSplitByStrand(opts.mVM.count("strand"));
	w.setGzip(opts.mVM.count("gzip"));

	OnlineBamMultiReportApp wigpe(tracer, &w);
	wigpe.processReads(input, output);

	tracer << "Total reads:\n";
	tracer << " On pos strand:\t" << w.getPosCnt() << "\n";
	tracer << " On neg strand:\t" << w.getNegCnt() << "\n";
	tracer << "\nWiggle file generation completed.\n\n";

}

std::string WigPE::normOutFileName(const std::string& fn) {
	std::string res(fn);
	if (boost::find_first(res, ".wig")) {
		erase_all(res, ".wig");
	}
	return res;
}

} /* namespace app */
