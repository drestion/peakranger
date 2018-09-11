/*
 * chip4cwig.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#include "chip4cwig.h"
#include "gzstream.h"

#include "region_profile/profilezoom.h"
#include "short_reads/readstools.h"
#include <boost/algorithm/string.hpp>
using namespace std;
using namespace boost;

chip_4c_wig::chip_4c_wig() {
    this->_readextlength = 200;
    this->_binlength = 100000000;
    this->_readlength = 25;
    this->_ws = 10000;
    this->_ov = 9000;
    this->_repeat = 5;
}

chip_4c_wig::~chip_4c_wig() {

}

void chip_4c_wig::smooth(wigs & _wigs, wigs & _r, uint32_t window, uint32_t ov) {
    profile_zoom _z;
    _z.smooth(window, ov, _wigs, _r);
}

void chip_4c_wig::export_smoothed_wig(Reads & reads, const char *file, bool gzip) {
    MARK_FUN( "chip_4c_wig::export_smoothed_wig");

    LOG_DEBUG1( "window size: " << _ws << "\n");LOG_DEBUG1( "overlap size:" << _ov << "\n");
    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart, preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart, npreadsend, nppreadsend, npppreadsend;
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);
    string sfile(file);

    boost::replace_last(sfile, ".wig", "");
    string pf(sfile + ".wig");

    ofstream pof(pf.c_str());

    rt_assert_msg(pof.good(), "file not good")

    _ct.print_msg(pof);

    print_wigfile_trackheader(pof, getPosRgb(), pf.c_str());
    wigs _wigs, _wigss;
    uint32_t _rp = _repeat;
    foreach(string chr, mergedchrs) {
        LOG_DEBUG3("In export_wiggle, start processing chr: "<<chr);
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

        pof << "variableStep chrom=" << chr << " span=1\n";

        _w._binned_wig_compiler(_binlength, _readlength, _readextlength, preadsstart, preadsend, npreadsstart,
                npreadsend, _wigss);

        while (_rp > 0) {

            smooth(_wigss, _wigs, _ws, _ov);

            _wigss.clear();

            std::copy(_wigs.begin(), _wigs.end(), std::back_inserter(_wigss));

            _wigs.clear();
            _rp--;
        }

        for (size_t i = 0; i < _wigss.size(); i++) {
            pof << _wigss[i].getP() << "\t" << _wigss[i].getS() << "\n";
        }

        _wigss.clear();
        _wigs.clear();
        _rp = _repeat;
    }
}

uint32_t chip_4c_wig::getBinlength() const {
    return _binlength;
}

uint32_t chip_4c_wig::getOv() const {
    return _ov;
}

uint32_t chip_4c_wig::getReadextlength() const {
    return _readextlength;
}

uint32_t chip_4c_wig::getReadlength() const {
    return _readlength;
}

wig_builder chip_4c_wig::getW() const {
    return _w;
}

uint32_t chip_4c_wig::getWs() const {
    return _ws;
}

void chip_4c_wig::setBinlength(uint32_t _binlength) {
    this->_binlength = _binlength;
}

void chip_4c_wig::setOv(uint32_t _ov) {
    this->_ov = _ov;
}

void chip_4c_wig::setReadextlength(uint32_t _readextlength) {
    this->_readextlength = _readextlength;
}

void chip_4c_wig::setReadlength(uint32_t _readlength) {
    this->_readlength = _readlength;
}

void chip_4c_wig::export_splitted_smoothed_wig(Reads& reads, const char* file, bool gzip) {

    MARK_FUN( "chip_4c_wig::export_splitted_smoothed_wig");

    LOG_DEBUG1( "window size: " << _ws << "\n");LOG_DEBUG1( "overlap size:" << _ov << "\n");
    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart, preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart, npreadsend, nppreadsend, npppreadsend;
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);
    string sfile(file);

    boost::replace_last(sfile, ".wig", "");

    foreach(string chr, mergedchrs) {
        string pf(sfile+"_"+chr+".wig");
        ofstream pof(pf.c_str());
        rt_assert_msg(pof.good(), "file not good")
        _ct.print_msg(pof);

        print_wigfile_trackheader(pof, getPosRgb(), pf.c_str());

        wigs _wigs, _wigss;
        uint32_t _rp = _repeat;
        LOG_DEBUG3("In export_wiggle, start processing chr: "<<chr);
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

        pof << "variableStep chrom=" << chr << " span=1\n";

        _w._binned_wig_compiler(_binlength, _readlength, _readextlength, preadsstart, preadsend, npreadsstart,
                npreadsend, _wigss);

        while (_rp > 0) {

            smooth(_wigss, _wigs, _ws, _ov);

            _wigss.clear();

            std::copy(_wigs.begin(), _wigs.end(), std::back_inserter(_wigss));

            _wigs.clear();
            _rp--;
        }

        for (size_t i = 0; i < _wigss.size(); i++) {
            pof << _wigss[i].getP() << "\t" << _wigss[i].getS() << "\n";
        }

        _wigss.clear();
        _wigs.clear();
        _rp = _repeat;
    }
}

void chip_4c_wig::setW(wig_builder _w) {
    this->_w = _w;
}

void chip_4c_wig::setWs(uint32_t _ws) {
    this->_ws = _ws;
}

uint32_t chip_4c_wig::getRepeat() const {
    return _repeat;
}

void chip_4c_wig::setRepeat(uint32_t _repeat) {
    this->_repeat = _repeat;
}

