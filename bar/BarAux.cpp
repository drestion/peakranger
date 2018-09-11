/*
 * BarAux.cpp
 *
 *  Created on: May 25, 2012
 *      Author: tania
 */

#include "BarAux.h"
#include "common/boost_header.h"
#include "common/stl_header.h"
#include "utils/stringutil.h"
using namespace boost::icl;
using namespace ranger::concepts;
using namespace std;
using namespace reads;

namespace ranger {
namespace bar {

void accHits(const std::vector<reads::Read>& hits, BarCounter& counter) {
	foreach(const Read& r, hits) {
		accHit(r, counter);
	}
}

void getCount(vector<RegionCount>& res, BarCounter& counter) {
	BarCounter::const_iterator it;
	for (it = counter.begin(); it != counter.end(); it++) {
		ICLInterval itv = (*it).first;
		int32_t cnt = (*it).second;
		res.push_back(RegionCount(RegionInt32(itv.lower(), itv.upper()), cnt));
	}
}

void accHit(const reads::Read& r, boost::icl::BarCounter& counter) {
	ICLInterval r1Val = ICLInterval::right_open(r.getStart(), r.getEnd());
	counter += make_pair(r1Val, 1);
}

void printBedGraphTrackline(std::ostream& os, bool addNewLine) {
	os << "track type=bedGraph autoScale=on";
	if (addNewLine)
		os << "\n";
}
void printBedGraphTrackline(ostream & pof, const char* _name,
		vector<uint32_t> col, bool isNewLine) {

	if (col.size() < 3) {
//            uint32_t ga[3] = { 225, 127, 0 };
		uint32_t ga[3] = { 50, 126, 184 };
		col = vector<uint32_t>(ga, ga + 3);
	}

	pof << "track type=bedGraph name=\"";
	pof << utils::ExtractFilename(string(_name));
	pof << "\" color=" << col[0] << "," << col[1] << "," << col[2];
	pof << " autoScale=on";
	if (isNewLine)
		pof << "\n";
}

}
}
