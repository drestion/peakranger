/*
 * BGHistStream.cpp
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#include "bar/BGHistStream.h"
#include "short_reads/Read.h"
#include "wiggle/gzstream.h"
#include "utils/Stamp.h"
using namespace boost::icl;
using namespace std;
using namespace utils;
namespace ranger {
namespace bar {

void outputAll(const std::map<string, BarCounter>& counter, ostream& os) {
    std::map<string, BarCounter>::const_iterator it;
    for (it = counter.begin(); it != counter.end(); ++it) {
        outputCounter(it->second, it->first, os);
    }
}

void outputCounter(const BarCounter& counter, const string& chr, ostream& os) {
    BarCounter::const_iterator it2;
    for (it2 = counter.begin(); it2 != counter.end(); it2++) {
        ICLInterval itv = (*it2).first;
        int32_t cnt = (*it2).second;
        os << chr << "\t" << itv.lower() << "\t" << itv.upper() << "\t" << cnt
                << "\n";
    }
}

void outputBar(const BarHist& hist, const string& chr, ostream& os) {
    const BarCounter& c = hist.counter();
    outputCounter(c, chr, os);
}

void combineAndOutput(const BarHist& pos, const BarHist& neg, const string& chr,
        ostream& os) {
    BarHist all(pos);
    all += neg;
    outputBar(all, chr, os);
}

void outputAll(const std::map<string, BarHist>& m, ostream& os) {
    std::map<string, BarHist>::const_iterator it;
    for (it = m.begin(); it != m.end(); ++it) {
        outputBar(it->second, it->first, os);
    }
}

void outputAllByChrWithHeader(std::map<string, BarHist>& counter,
        const string& file_prefix, const string& file_postfix, bool gzip) {
    std::map<string, BarHist>::iterator it;
    for (it = counter.begin(); it != counter.end(); ++it) {
        string file(file_prefix + "_");
        file += it->first;
        file += file_postfix;

        if (gzip) {
            file += ".gz";
            ogzstream os(file.c_str());
            Stamp::citationAndDate(os);
            printBedGraphTrackline(os, file.c_str());
            outputBar(it->second, it->first, os);
        } else {
            ofstream os(file.c_str());
            Stamp::citationAndDate(os);
            printBedGraphTrackline(os, file.c_str());
            outputBar(it->second, it->first, os);
        }
    }
}

void outputAllWithHeader(const std::map<string, BarHist>& counter,
        const string& file_prefix, bool gzip) {
    if (gzip) {
        string f(file_prefix);
        f += ".gz";
        ogzstream os(f.c_str());
        Stamp::citationAndDate(os);
        printBedGraphTrackline(os, f.c_str());
        outputAll(counter, os);
    } else {
        ofstream os(file_prefix.c_str());
        Stamp::citationAndDate(os);
        printBedGraphTrackline(os, file_prefix.c_str());
        outputAll(counter, os);
    }
}

void combineAndOutput(const BarHist& pos, const BarHist& neg, const string& chr,
        const string& file, bool gzip) {
    BarHist all(pos);
    all += neg;
    if (gzip) {
        string f(file);
        f += ".gz";
        ogzstream os(f.c_str());
        Stamp::citationAndDate(os);
        outputBar(all, chr, os);
    } else {
        ofstream os(file.c_str());
        Stamp::citationAndDate(os);
        outputBar(all, chr, os);
    }

}


} /* namespace bar */
} /* namespace ranger */
