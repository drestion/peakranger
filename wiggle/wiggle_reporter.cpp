/*
 * wiggle_reporter.cpp
 *
 *  Created on: Jun 9, 2011
 *      Author: xin
 */

#include "wiggle_reporter.h"
#include "region_profile/RegionProfile.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include "utils/stringutil.h"

#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <ostream>
using namespace std;
using namespace boost;
using namespace utils;
namespace {
bool wiggle_sgr_sort_comparator(pair<uint32_t, double> p1
                                ,
                                pair<uint32_t, double> p2) {
    return p1.first < p2.first;
}
}
void wiggle_reporter::split_export_wiggle(Reads & reads,
                                          ostream & os) {
    throw RangerException("Sorry, wiggle_reporter::split_export_wiggle(Reads "
                          "& reads,ostream & os) not implemented yet.");
}

void wiggle_reporter::split_export_wiggle(Reads & reads,
                                          const char *file) {
    string sfile(file);
    size_t ind = sfile.find_last_of(".wig");
    if (ind != string::npos) {
        //remove .wig
        sfile = sfile.substr(0,
                             ind);
    }

    assert_eq(reads.pos_reads.chrs().size(),
              reads.neg_reads.chrs().size())

    foreach(string chr, reads.pos_reads.chrs())
    {

        setWiggleName(sfile + chr + ".wig");
        string newfile(sfile + chr + ".wig");
        ofstream ofs(newfile.c_str());
        assert(ofs.is_open());
        if (!(ofs.is_open())) throw FileNotGood(newfile.c_str());

        vector<uint32_t> preads(reads.pos_reads.begin_of(chr),
                                reads.pos_reads.end_of(chr));
        vector<uint32_t> nreads(reads.neg_reads.begin_of(chr),
                                reads.neg_reads.end_of(chr));

        export_wiggle(preads,
                      nreads,
                      chr,
                      ofs);
        ofs.close();
    }
}

void wiggle_reporter::export_wiggle(vector<uint32_t> & preads,
                                    vector<uint32_t> & nreads,
                                    string chr,
                                    ostream & os) {

    os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
    << "visibility=dense " << "color=" << _colorRGB[0] << "," << _colorRGB[1]
    << "," << _colorRGB[2] << " " << "altColor=" << _colorRGB[0] << ","
    << _colorRGB[1] << "," << _colorRGB[2] << " " << "priority=" << _priority
    << "\n";

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
    preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
    npreadsend, nppreadsend, npppreadsend;
    assert_neq(chr,
               "") LOG_DEBUG1("In export_wiggle, start processing chr: "<<chr);
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
    uint32_t binstart = 0;
    uint32_t binend = 0;

    os << "variableStep chrom=" << chr << " span=1\n";

    assert_gt(_binlength,
              1)
    while (noofbins-- > 0) {

        LOG_DEBUG2("start Bin: "<<binind);
        binstart = _binlength * binind + 1;
        binend = _binlength * (binind + 1);
        binind++;

        ppreadsstart = lower_bound(preadsstart,
                                   preadsend,
                                   binstart);
        ppreadsend = upper_bound(ppreadsstart,
                                 preadsend,
                                 binend);
        nppreadsstart = lower_bound(npreadsstart,
                                    npreadsend,
                                    binstart);
        nppreadsend = upper_bound(nppreadsstart,
                                  npreadsend,
                                  binend);

        this->_process(binstart,
                       binend,
                       _readlength,
                       _readextlength,
                       ppreadsstart,
                       ppreadsend,
                       nppreadsstart,
                       nppreadsend,
                       os);
        LOG_DEBUG2("Bin processed: "<<binind);

    }LOG_DEBUG1("In export_wiggle, finished processing chr: "<<chr);

}

void wiggle_reporter::print_wig(data_type& reads_count,
                                ostream& os) {
    for (size_t i = 0; i < reads_count.size(); i++) {
        os << reads_count[i].first << "\t" << reads_count[i].second << "\n";
    }
}

void wiggle_reporter::_process(uint32_t start,
                               uint32_t end,
                               uint32_t readlength,
                               uint32_t readextlength,
                               vector<uint32_t>::iterator readsStart,
                               vector<uint32_t>::iterator readsEnd,
                               vector<uint32_t>::iterator nreadsStart,
                               vector<uint32_t>::iterator nreadsEnd,
                               ostream& os) {
    assert_gt(end,
              2)

    assert_gt(end - 2,
              start)

    assert_gt(end-start+1,
              readlength)

    assert_gt(end-start+1,
              readextlength)

    LOG_DEBUG1("wiggle_reporter::_process");

    LOG_DEBUG1("readlength:"<<readlength);

    LOG_DEBUG1("readextlength:"<<readextlength);
    uint32_t a;
    uint32_t b;
    uint32_t arrayStart;
    uint32_t arrayEnd;
    uint32_t read;
    typedef vector<pair<uint32_t, double> > reads_count_t;
    reads_count_t reads_count;

    bool inRange = false;
    LOG_DEBUG1("start mapping pos reads");

    LOG_DEBUG1("Total pos reads: "<<readsEnd-readsStart);
    while (readsStart != readsEnd) {

        read = *readsStart;
        readsStart++;
        a = read;
        b = read + readextlength;
        arrayStart = 0;
        arrayEnd = 0;

        if (a < end && b > start) {
            inRange = true;
            /*
             *     |-------|
             *   |---|
             */
            if (a <= start && b <= end) {

                reads_count.push_back(pair<uint32_t, double>(start,
                                                             1));
                reads_count.push_back(pair<uint32_t, double>(b,
                                                             -1));
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {

                reads_count.push_back(pair<uint32_t, double>(a,
                                                             1));
                reads_count.push_back(pair<uint32_t, double>(b,
                                                             -1));
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {

                reads_count.push_back(pair<uint32_t, double>(a,
                                                             1));
                reads_count.push_back(pair<uint32_t, double>(end,
                                                             1));
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
        if (read + readlength < readextlength) continue;
        a = read + readlength - readextlength;
        b = read;
        arrayStart = 0;
        arrayEnd = 0;

        if (a < end && b > start) {
            inRange = true;
            /*
             *     |-------|
             *   |---|
             */
            if (a <= start && b <= end) {

                reads_count.push_back(pair<uint32_t, double>(start,
                                                             1));
                reads_count.push_back(pair<uint32_t, double>(b,
                                                             -1));
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {

                reads_count.push_back(pair<uint32_t, double>(a,
                                                             1));
                reads_count.push_back(pair<uint32_t, double>(b,
                                                             -1));
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {

                reads_count.push_back(pair<uint32_t, double>(a,
                                                             1));
                reads_count.push_back(pair<uint32_t, double>(end,
                                                             1));
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
    //sort based on their locations;
    sort(reads_count.begin(),
         reads_count.end(),
         wiggle_sgr_sort_comparator);

    uint32_t pos;
    uint32_t ppos = (reads_count.begin())->first;
    double count = (reads_count.begin())->second;
    LOG_DEBUG1("start building sgr");LOG_DEBUG1("initial count value:"<<count);
    for (size_t i = 1; i < reads_count.size(); i++) {
        pos = reads_count[i].first;
//        cout << "ps score:" << (double) (reads_count[i].second) << "\tps pos:" << pos << endl;
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

    LOG_DEBUG1("QUIT: wiggle_reporter::_process");
}

void wiggle_reporter::export_wiggle(Reads& reads,
                                    ostream & os) {

    os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
    << "visibility=dense autoScale=off " << "viewLimits=" << _viewLimitDown
    << ":" << _viewLimitUp << " " << "color=" << _colorRGB[0] << ","
    << _colorRGB[1] << "," << _colorRGB[2] << " " << "altColor=" << _colorRGB[0]
    << "," << _colorRGB[1] << "," << _colorRGB[2] << " " << "priority="
    << _priority << "\n";

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
    preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
    npreadsend, nppreadsend, npppreadsend;

    map<string, vector<uint32_t> >::iterator rit;
    rit = reads.pos_reads.begin();
    string chr;
    assert_eq(reads.pos_reads.chrs().size(),
              reads.neg_reads.chrs().size())
    for (; rit != reads.pos_reads.end(); rit++) {
        chr = rit->first;

        LOG_DEBUG1("In export_wiggle, start processing chr: "<<chr);
        preadsstart = reads.pos_reads.begin_of(chr);
        preadsend = reads.pos_reads.end_of(chr);
        npreadsstart = reads.neg_reads.begin_of(chr);
        npreadsend = reads.neg_reads.end_of(chr);

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
                              npreadsend == npreadsstart ? 0 :
                              (*(npreadsend - 1));

        uint32_t apchrlength =
                               pchrlength > nchrlength ? pchrlength :
                               nchrlength;
        uint32_t noofbins = apchrlength ? 1 + (apchrlength / _binlength) : 0;
        uint32_t binind = 0;
        uint32_t binstart, binend;

        os << "variableStep chrom=" << chr << " span=1\n";
        assert_gt(_binlength,
                  1);
        while (noofbins-- > 0) {
            LOG_DEBUG2("start Bin: "<<binind);
            binstart = _binlength * binind + 1;
            binend = _binlength * (binind + 1);
            binind++;

            ppreadsstart = lower_bound(preadsstart,
                                       preadsend,
                                       binstart);
            ppreadsend = upper_bound(ppreadsstart,
                                     preadsend,
                                     binend);
            nppreadsstart = lower_bound(npreadsstart,
                                        npreadsend,
                                        binstart);
            nppreadsend = upper_bound(nppreadsstart,
                                      npreadsend,
                                      binend);

            _process(binstart,
                     binend,
                     reads.getReadlength(),
                     _readextlength,
                     ppreadsstart,
                     ppreadsend,
                     nppreadsstart,
                     nppreadsend,
                     os);
            LOG_DEBUG2("Bin processed: "<<binind);

        }LOG_DEBUG1("In export_wiggle, finished processing chr: "<<chr);
    }
}

vector<uint32_t> wiggle_reporter::getColorRGB() const {
    return _colorRGB;
}

uint32_t wiggle_reporter::getPriority() const {
    return _priority;
}

uint32_t wiggle_reporter::getViewLimitDown() const {
    return _viewLimitDown;
}

uint32_t wiggle_reporter::getViewLimitUp() const {
    return _viewLimitUp;
}

string wiggle_reporter::getWiggleName() const {
    return _wiggleName;
}

void wiggle_reporter::setColorRGB(vector<uint32_t> _colorRGB) {
    assert_gt(_colorRGB.size(),
              3)assert_lt(_colorRGB.at(0),
                          256)assert_lt(_colorRGB.at(1),
                                        256)assert_lt(_colorRGB.at(2),
                                                      256)
    this->_colorRGB = _colorRGB;

}

void wiggle_reporter::setPriority(uint32_t _priority) {
    this->_priority = _priority;
}

void wiggle_reporter::setViewLimitDown(uint32_t _viewLimitDown) {
    this->_viewLimitDown = _viewLimitDown;
}

void wiggle_reporter::setViewLimitUp(uint32_t _viewLimitUp) {
    this->_viewLimitUp = _viewLimitUp;
}

void wiggle_reporter::export_wiggle_gzip(Reads & reads,
                                         const char *file) {
}

void wiggle_reporter::split_export_wiggle_gzip(Reads & reads,
                                               const char *file) {
}

void wiggle_reporter::setWiggleName(string _wiggleName) {
    assert_neq(_wiggleName,
               "");
    string sfile(_wiggleName);
    string dir, ext, file;
    stringutil::get_dir_file(_wiggleName,
                             dir,
                             file,
                             ext);
    string filename = "PeakRanger_" + file;
    size_t ind = filename.find_last_of(".wig");
    if (ind == string::npos) {
        //add .wig
        filename += ".wig";
    }
    this->_wiggleName = filename;
}

void wiggle_reporter::export_wiggle(Reads & reads,
                                    const char *file) {

    string sfile(file);
    size_t ind = sfile.find_last_of(".wig");
    if (ind == string::npos) {
        //add .wig
        sfile += ".wig";
    }
    setWiggleName(string(file));

    ofstream ofs(sfile.c_str());
    if (ofs.is_open()) {
        export_wiggle(reads,
                      ofs);
        ofs.close();
    } else {
        throw FileNotGood(sfile.c_str());
    }
}

void wiggle_reporter::use_default_setting() {
    //    _wiggleName = "Positive reads count.PeakRanger";
    uint32_t ga[3] = { 50, 126, 184 };
    _colorRGB = vector < uint32_t > (ga,
    ga + 3);
    _priority = 20;
    _viewLimitDown = 0;
    _viewLimitUp = 1024;
    _binlength = 100000000;
    _readextlength = 100;
    _readlength = 25;
}

uint32_t wiggle_reporter::getBinlength() const {
    return _binlength;
}

void wiggle_reporter::setBinlength(uint32_t _binlength) {
    this->_binlength = _binlength;
}

uint32_t wiggle_reporter::getReadextlength() const {
    return _readextlength;
}

void wiggle_reporter::setReadextlength(uint32_t _readextlength) {
    this->_readextlength = _readextlength;
}

uint32_t wiggle_reporter::getReadlength() const {
    return _readlength;
}

void wiggle_reporter::setReadlength(uint32_t _readlength) {
    this->_readlength = _readlength;
}

