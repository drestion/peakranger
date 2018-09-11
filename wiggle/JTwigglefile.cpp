/*
 * JTwigglefile.cpp
 *
 *  Created on: Oct 4, 2011
 *      Author: xfeng
 */

#include "JTwigglefile.h"
#include "utils/logger.h"
#include "short_reads/readstools.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include "utils/stringutil.h"
#include "utils/Stamp.h"
#include "gzstream.h"

#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
using namespace std;
using namespace boost;

namespace {
bool wiggle_sgr_sort_comparator(pair<uint32_t, double> p1,
		pair<uint32_t, double> p2) {
	return p1.first < p2.first;
}
}

JT_wiggle_file::JT_wiggle_file() {

	this->setViewLimitUp((uint32_t) 512);
}

JT_wiggle_file::~JT_wiggle_file() {

}

void JT_wiggle_file::_process(uint32_t start, uint32_t end, uint32_t readlength,
		uint32_t readextlength, vector<uint32_t>::iterator readsStart,
		vector<uint32_t>::iterator readsEnd,
		vector<uint32_t>::iterator nreadsStart,
		vector<uint32_t>::iterator nreadsEnd, ostream& os) {
//    cout << "in JT_process" << endl;
	assert_gt(end, 2)
	assert_gt(end - 2, start)
	assert_gt(end - start + 1, readlength)
	assert_gt(end - start + 1, readextlength) LOG_DEBUG1("JT_JT_wiggle_file::_process");LOG_DEBUG1("readlength:"<<readlength);LOG_DEBUG1("readextlength:"<<readextlength);
	uint32_t a;
	uint32_t b;
	uint32_t arrayStart;
	uint32_t arrayEnd;
	uint32_t read;
	typedef vector<pair<uint32_t, double> > reads_count_t;
	reads_count_t reads_count;

	bool inRange = false;
	LOG_DEBUG1("start mapping pos reads");LOG_DEBUG1("Total pos reads: "<<readsEnd-readsStart);
	while (readsStart != readsEnd) {

		read = *readsStart;
		readsStart++;
		a = read;
		b = read + readextlength;
//        cout <<"read:"<<read<< " a:"<<a<<" b:"<<b<<endl;
		arrayStart = 0;
		arrayEnd = 0;

		if (a < end && b > start) {
			inRange = true;
			/*
			 *     |-------|
			 *   |---|
			 */
			if (a <= start && b <= end) {

				reads_count.push_back(pair<uint32_t, double>(start, 1));
				reads_count.push_back(pair<uint32_t, double>(b, -1));
			}
			/*
			 *     |-------|
			 *       |---|
			 */
			if (a >= start && b <= end) {

				reads_count.push_back(pair<uint32_t, double>(a, 1));
				reads_count.push_back(pair<uint32_t, double>(b, -1));
			}
			/*
			 *     |-------|
			 *            |---|
			 */
			if (a < end && b > end) {

				reads_count.push_back(pair<uint32_t, double>(a, 1));
				reads_count.push_back(pair<uint32_t, double>(end, 1));
			}
			/*
			 *     |-------|
			 *   |-----------|
			 */
			if (a < start && b > end) {

			}

		} else if (inRange) {
			break;
		}
	}

	LOG_DEBUG1("start mapping neg reads");LOG_DEBUG1("Total neg reads: "<<nreadsEnd-nreadsStart);
	while (nreadsStart != nreadsEnd) {

		read = *nreadsStart;
		nreadsStart++;
		if (read + readlength < readextlength)
			continue;
		// a = read + readlength - readextlength;
//        b = read;
		//THIS IS WHERE JT style kicks in

		a = read + readlength - readextlength;
		b = read + readlength;
//        cout <<"wig negread:"<<read<< " a:"<<a<<" b:"<<b<<endl;
		arrayStart = 0;
		arrayEnd = 0;

		if (a < end && b > start) {
			inRange = true;
			/*
			 *     |-------|
			 *   |---|
			 */
			if (a <= start && b <= end) {

				reads_count.push_back(pair<uint32_t, double>(start, 1));
				reads_count.push_back(pair<uint32_t, double>(b, -1));
			}
			/*
			 *     |-------|
			 *       |---|
			 */
			if (a >= start && b <= end) {

				reads_count.push_back(pair<uint32_t, double>(a, 1));
				reads_count.push_back(pair<uint32_t, double>(b, -1));
			}
			/*
			 *     |-------|
			 *            |---|
			 */
			if (a < end && b > end) {

				reads_count.push_back(pair<uint32_t, double>(a, 1));
				reads_count.push_back(pair<uint32_t, double>(end, 1));
			}
			/*
			 *     |-------|
			 *   |-----------|
			 */
			if (a < start && b > end) {

			}

		} else if (inRange) {
			break;
		}
	}
	if (reads_count.size() == 0) {
		return;
	}
	if (reads_count.size())
		//sort based on their locations;
		sort(reads_count.begin(), reads_count.end(),
				wiggle_sgr_sort_comparator);

	uint32_t pos;
	uint32_t ppos = (reads_count.begin())->first;
	double count = (reads_count.begin())->second;
	LOG_DEBUG1("start building sgr");LOG_DEBUG1("initial count value:"<<count);

	for (size_t i = 1; i < reads_count.size(); i++) {
		pos = reads_count[i].first;
//        cout <<"ps score:"<<(double)(reads_count[i].second)<<"\tps pos:"<<pos<<endl;
		if (pos == ppos) {
			count += (double) (reads_count[i].second);
			continue;
		}
		if (ppos > 0 && count > 0) {
			//This info is viewable in the result wig file
			os << ppos << "\t" << count << "\n";
		}
		count += (double) (reads_count[i].second);
		ppos = pos;
	}

	if (ppos > 0 && count > 0) {
		os << ppos << "\t" << count << "\n";
	}

	LOG_DEBUG1("QUIT: JT_JT_wiggle_file::_process");
}

void JT_wiggle_file::split_export_wiggle(Reads & reads, ostream & os) {
	throw RangerException("Sorry, JT_wiggle_file::split_export_wiggle(Reads "
			"& reads,ostream & os) not implemented yet.");
}

void JT_wiggle_file::export_wiggle(vector<uint32_t> & preads,
		vector<uint32_t> & nreads, string chr, ostream & os) {
//    cout << "In jT_sp_export_wiggle" << endl;
	assert_neq(chr, "")
	os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
			<< "visibility=dense " << "color=" << _colorRGB[0] << ","
			<< _colorRGB[1] << "," << _colorRGB[2] << " " << "altColor="
			<< _colorRGB[0] << "," << _colorRGB[1] << "," << _colorRGB[2] << " "
			<< "priority=" << _priority << "\n";

	vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
			preadsend, ppreadsend, pppreadsend;

	vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
			npreadsend, nppreadsend, npppreadsend;

	LOG_DEBUG1("In export_wiggle, start processing chr: "<<chr);
	preadsstart = preads.begin();
	preadsend = preads.end();
	npreadsstart = nreads.begin();
	npreadsend = nreads.end();

	/*
	 * in case some chrs only contain pos or neg reads
	 * do not use
	 * pchrlength = (*(preadsend-1));
	 */
	//ignore the strand if it contains less than 2 read
	//todo: this can not rule out the case preadsend = 0x00 if this function
	// is not called after reads correction.
	uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
	uint32_t nchrlength = npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));

	uint32_t apchrlength = pchrlength > nchrlength ? pchrlength : nchrlength;
	uint32_t noofbins = apchrlength ? 1 + (apchrlength / _binlength) : 0;
	uint32_t binind = 0;
	uint32_t binstart, binend;

	os << "variableStep chrom=" << chr << " span=1\n";
	assert_gt(_binlength, 1)
	while (noofbins-- > 0) {
		LOG_DEBUG2("start Bin: "<<binind);
		binstart = _binlength * binind + 1;
		binend = _binlength * (binind + 1);
		binind++;

		ppreadsstart = lower_bound(preadsstart, preadsend, binstart);
		ppreadsend = upper_bound(ppreadsstart, preadsend, binend);
		nppreadsstart = lower_bound(npreadsstart, npreadsend, binstart);
		nppreadsend = upper_bound(nppreadsstart, npreadsend, binend);

		this->_process(binstart, binend, _readlength, _readextlength,
				ppreadsstart, ppreadsend, nppreadsstart, nppreadsend, os);
		LOG_DEBUG2("Bin processed: "<<binind);

	}LOG_DEBUG1("In export_wiggle, finished processing chr: "<<chr);

}
void JT_wiggle_file::export_wiggle(Reads & reads, ostream & os) {

	os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
			<< "visibility=dense " << "color=" << _colorRGB[0] << ","
			<< _colorRGB[1] << "," << _colorRGB[2] << " " << "altColor="
			<< _colorRGB[0] << "," << _colorRGB[1] << "," << _colorRGB[2] << " "
			<< "priority=" << _priority << "\n";

	vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
			preadsend, ppreadsend, pppreadsend;

	vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
			npreadsend, nppreadsend, npppreadsend;
	vector<string> mergedchrs;
	reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);

	foreach(string chr, mergedchrs) {
		LOG_DEBUG1("In export_wiggle, start processing chr: "<<chr);
		if (reads.pos_reads.hasReadsOnChr(chr)) {
			preadsstart = reads.pos_reads.begin_of(chr);
			preadsend = reads.pos_reads.end_of(chr);
		} else {
			preadsstart = preadsend;
		}
		if (reads.neg_reads.hasReadsOnChr(chr)) {
			npreadsstart = reads.neg_reads.begin_of(chr);
			npreadsend = reads.neg_reads.end_of(chr);
		} else {
			npreadsend = npreadsstart;
		}

		/*
		 * in case some chrs only contain pos or neg reads
		 * do not use
		 * pchrlength = (*(preadsend-1));
		 */
		//ignore the strand if it contains less than 2 read
		//todo: this can not rule out the case preadsend = 0x00 if this function
		// is not called after reads correction.
		uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
		uint32_t nchrlength =
				npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));

		uint32_t apchrlength =
				pchrlength > nchrlength ? pchrlength : nchrlength;
		uint32_t noofbins = apchrlength ? 1 + (apchrlength / _binlength) : 0;
		uint32_t binind = 0;
		uint32_t binstart, binend;
		LOG_DEBUG2("pchrlength: "<<pchrlength<<" nchrlength:"<<nchrlength);
		os << "variableStep chrom=" << chr << " span=1\n";
		assert_gt(_binlength, 1)
		while (noofbins-- > 0) {
			LOG_DEBUG2("start Bin: "<<binind);
			binstart = _binlength * binind + 1;
			binend = _binlength * (binind + 1);
			LOG_DEBUG2("binstart: "<<binstart<<" binend:"<<binend);
			binind++;

			if (preadsstart != preadsend) {
				ppreadsstart = lower_bound(preadsstart, preadsend, binstart);
				ppreadsend = upper_bound(ppreadsstart, preadsend, binend);
			} else {
				ppreadsstart = ppreadsend;
			}
			if (npreadsstart != npreadsend) {
				nppreadsstart = lower_bound(npreadsstart, npreadsend, binstart);
				nppreadsend = upper_bound(nppreadsstart, npreadsend, binend);
			} else {
				nppreadsstart = nppreadsend;
			}
			_process(binstart, binend, reads.getReadlength(), _readextlength,
					ppreadsstart, ppreadsend, nppreadsstart, nppreadsend, os);
			LOG_DEBUG2("Bin processed: "<<binind);

		}LOG_DEBUG1("In export_wiggle, finished processing chr: "<<chr);
	}
}
void JT_wiggle_file::split_export_wiggle(Reads & reads, const char *file) {
//    cout << "In master jT_sp_export_wiggle" << endl;
	LOG_DEBUG1("JT_wiggle_file::split_export_wiggle");
	rt_assert(file)
	string sfile(file);
	size_t ind = sfile.find_last_of(".wig");
	if (ind != string::npos) {
		//remove .wig
		sfile = sfile.substr(0, ind - 3);
	}
	vector<string> mergedchrs;
	reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);
#ifdef USE_LOGGING
	foreach(string chr, reads.pos_reads.chrs()) {
		LOG_DEBUG1("POS chr:"<<chr);
	}
	foreach(string chr, reads.neg_reads.chrs()) {
		LOG_DEBUG1("Neg chr:"<<chr);
	}
	foreach(string chr, mergedchrs) {
		LOG_DEBUG1("Merged chr:"<<chr);
	}
#endif
	foreach(string chr, mergedchrs) {
		setWiggleName(sfile + "_" + chr + ".wig");

		string newfile(sfile +"_" + chr + ".wig");
		ofstream ofs(newfile.c_str());
		utils::Stamp::citationAndDate(ofs);
		rt_assert_msg(ofs.is_open(), "file not good")
		vector<uint32_t> preads;

		vector<uint32_t> nreads;
		if (reads.pos_reads.hasReadsOnChr(chr)) {
			preads.insert(preads.begin(), reads.pos_reads.begin_of(chr),
					reads.pos_reads.end_of(chr));
		}

		if (reads.neg_reads.hasReadsOnChr(chr)) {
			nreads.insert(nreads.begin(), reads.neg_reads.begin_of(chr),
					reads.neg_reads.end_of(chr));
		}

		LOG_DEBUG1("Generating: "<<newfile);
		export_wiggle(preads, nreads, chr, ofs);
		ofs.close();
	}
}

void JT_wiggle_file::export_wiggle(Reads & reads, const char *file) {

	string sfile(file);
	size_t ind = sfile.find_last_of(".wig");
	if (ind == string::npos) {
		//add .wig
		sfile += ".wig";
	}
	setWiggleName(string(file));

	ofstream ofs(sfile.c_str());
	utils::Stamp::citationAndDate(ofs);
	rt_assert_msg(ofs.is_open(), "file not good")

	export_wiggle(reads, ofs);

	ofs.close();
}

void JT_wiggle_file::export_wiggle_gzip(Reads & reads, const char *file) {
	string sfile(file);
	size_t ind = sfile.find_last_of(".wig");
	if (ind != string::npos) {
		//remove .wig
		sfile = sfile.substr(0, ind - 3);
	}
	sfile += ".wig.gz";
	setWiggleName(sfile);

	ogzstream ofs(sfile.c_str());
	utils::Stamp::citationAndDate(ofs);
	rt_assert_msg(ofs.is_open(), "file not good")

	export_wiggle(reads, ofs);

	ofs.close();
}

void JT_wiggle_file::split_export_wiggle_gzip(Reads & reads, const char *file) {
	//    cout << "In master jT_sp_export_wiggle" << endl;
	LOG_DEBUG1("JT_wiggle_file::split_export_wiggle");
	rt_assert(file)
	string sfile(file);
	size_t ind = sfile.find_last_of(".wig");
	if (ind != string::npos) {
		//remove .wig
		sfile = sfile.substr(0, ind - 3);
	}
	vector<string> mergedchrs;
	reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);
#ifdef USE_LOGGING
	foreach(string chr, reads.pos_reads.chrs()) {
		LOG_DEBUG1("POS chr:"<<chr);
	}
	foreach(string chr, reads.neg_reads.chrs()) {
		LOG_DEBUG1("Neg chr:"<<chr);
	}
	foreach(string chr, mergedchrs) {
		LOG_DEBUG1("Merged chr:"<<chr);
	}
#endif
	foreach(string chr, mergedchrs) {
		string newfile(sfile + "_"+chr + ".wig.gz");
		setWiggleName(newfile);

		ogzstream ofs(newfile.c_str());
		utils::Stamp::citationAndDate(ofs);
		rt_assert_msg(ofs.is_open(), "file not good")
		vector<uint32_t> preads;

		vector<uint32_t> nreads;
		if (reads.pos_reads.hasReadsOnChr(chr)) {
			preads.insert(preads.begin(), reads.pos_reads.begin_of(chr),
					reads.pos_reads.end_of(chr));
		}

		if (reads.neg_reads.hasReadsOnChr(chr)) {
			nreads.insert(nreads.begin(), reads.neg_reads.begin_of(chr),
					reads.neg_reads.end_of(chr));
		}

		LOG_DEBUG1("Generating: "<<newfile);
		export_wiggle(preads, nreads, chr, ofs);
		ofs.close();
	}
}

