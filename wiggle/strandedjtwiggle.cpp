/*
 * strandedjtwiggle.cpp
 *
 *  Created on: Jan 11, 2012
 *      Author: xfeng
 */

#include "strandedjtwiggle.h"
#include "gzstream.h"
#include "utils/logger.h"
#include "short_reads/readstools.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include "utils/stringutil.h"
#include "utils/Stamp.h"
#include "region_profile/profilezoom.h"

#include <fstream>
#include <iostream>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

using namespace std;

namespace {
typedef pair<uint32_t, double> element_type;
bool wiggle_sgr_sort_comparator(element_type p1, element_type p2) {
    return p1.first < p2.first;
}
}
stranded_jtwiggle::stranded_jtwiggle() {
    uint32_t ga[3] = { 228, 26, 28 };
    _ncolorRGB = vector<uint32_t>(ga, ga + 3);

}

stranded_jtwiggle::~stranded_jtwiggle() {

}

void stranded_jtwiggle::export_wiggle(Reads & reads, ostream & os) {
    cout << "stranded_jtwiggle::export_wiggle is not conceptually right\n";
    exit(0);
}

void stranded_jtwiggle::export_wiggle(vector<uint32_t> & preads,
        vector<uint32_t> & nreads, string chr, ostream & os) {
    cout << "stranded_jtwiggle::export_wiggle is not conceptually right\n";
    exit(0);
}

void stranded_jtwiggle::export_wiggle(Reads & reads, const char *file) {
    string sfile(file);
    size_t ind = sfile.find_last_of(".wig");
    if (ind != string::npos) {
        //remove .wig
        sfile = sfile.substr(0, ind - 3);
    }

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
            preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
            npreadsend, nppreadsend, npppreadsend;
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);
    string pf(sfile + "_Pos.wig");
    string nf(sfile + "_Neg.wig");

    ofstream pof(pf.c_str());
    ofstream nof(nf.c_str());
    utils::Stamp::citationAndDate(pof);
    utils::Stamp::citationAndDate(nof);
    rt_assert_msg(pof.is_open(), "file not good")
    rt_assert_msg(nof.is_open(), "file not good")
    pof << "track type=wiggle_0 name=\"" << pf << "\" " << "visibility=dense "
            << "color=" << _colorRGB[0] << "," << _colorRGB[1] << ","
            << _colorRGB[2] << " " << "altColor=" << _colorRGB[0] << ","
            << _colorRGB[1] << "," << _colorRGB[2] << " " << "priority="
            << _priority << "\n";
    nof << "track type=wiggle_0 name=\"" << nf << "\" " << "visibility=dense "
            << "color=" << _ncolorRGB[0] << "," << _ncolorRGB[1] << ","
            << _ncolorRGB[2] << " " << "altColor=" << _ncolorRGB[0] << ","
            << _ncolorRGB[1] << "," << _ncolorRGB[2] << " " << "priority="
            << _priority << "\n";
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

        uint32_t noofbins = pchrlength ? 1 + (pchrlength / _binlength) : 0;
        uint32_t binind = 0;
        uint32_t binstart = 0, binend = 0;

        uint32_t nnoofbins = nchrlength ? 1 + (nchrlength / _binlength) : 0;
        uint32_t nbinind = 0;
        uint32_t nbinstart = 0, nbinend = 0;
        pof << "variableStep chrom=" << chr << " span=1\n";
        nof << "variableStep chrom=" << chr << " span=1\n";
        assert_gt(_binlength, 1)
        while (noofbins-- > 0) {
            LOG_DEBUG2("start Bin: "<<binind);
            binstart = _binlength * binind + 1;
            binend = _binlength * (binind + 1);
            binind++;

            if (preadsstart != preadsend) {
                ppreadsstart = lower_bound(preadsstart, preadsend, binstart);
                ppreadsend = upper_bound(ppreadsstart, preadsend, binend);
            } else {
                ppreadsstart = ppreadsend;
            }

            nppreadsstart = nppreadsend;

            JT_wiggle_file::_process(binstart, binend, reads.getReadlength(),
                    _readextlength, ppreadsstart, ppreadsend, nppreadsstart,
                    nppreadsend, pof);
            LOG_DEBUG2("pos Bin processed: "<<binind);

        }

        LOG_DEBUG1("In export_wiggle, finished processing chr pos: "<<chr);
        while (nnoofbins-- > 0) {
            LOG_DEBUG2("start Bin: "<<nbinind);
            nbinstart = _binlength * nbinind + 1;
            nbinend = _binlength * (nbinind + 1);
            nbinind++;

            ppreadsstart = ppreadsend;

            if (npreadsstart != npreadsend) {
                nppreadsstart = lower_bound(npreadsstart, npreadsend,
                        nbinstart);
                nppreadsend = upper_bound(nppreadsstart, npreadsend, nbinend);
            } else {
                nppreadsstart = nppreadsend;
            }
            _process_neg(nbinstart, nbinend, reads.getReadlength(),
                    _readextlength, nppreadsstart, nppreadsend, nof);
            LOG_DEBUG2("neg Bin processed: "<<nbinind);

        }

        LOG_DEBUG1("In export_wiggle, finished processing chr neg: "<<chr);
    }
    pof.close();
    nof.close();
}

void stranded_jtwiggle::split_export_wiggle(Reads & reads, ostream & os) {
    cout << "stranded_jtwiggle::split_export_wiggle not implemented yet.\n ";
    exit(0);
}

void stranded_jtwiggle::split_export_wiggle(Reads & reads, const char *file) {

    LOG_DEBUG1("stranded_jtwiggle::split_export_wiggle");
    rt_assert(file)
    string sfile(file);
    size_t ind = sfile.find_last_of(".wig");
    if (ind != string::npos) {
        //remove .wig
        sfile = sfile.substr(0, ind - 3);
    }
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);

    foreach(string chr, mergedchrs) {

        if (reads.pos_reads.hasReadsOnChr(chr)) {
            vector<uint32_t> preads;

            string newfile(sfile + "_" + chr + "_Pos.wig");
            preads.insert(preads.begin(), reads.pos_reads.begin_of(chr),
                    reads.pos_reads.end_of(chr));

            setWiggleName(newfile);

            ofstream ofs(newfile.c_str());
            utils::Stamp::citationAndDate(ofs);
            rt_assert_msg(ofs.is_open(), "file not good") LOG_DEBUG1("Generating: "<<newfile);
            export_wiggle_pos(preads, chr, ofs);
            ofs.close();
        }

        if (reads.neg_reads.hasReadsOnChr(chr)) {

            vector<uint32_t> nreads;
            string newfile(sfile + "_" + chr + "_Neg.wig");
            nreads.insert(nreads.begin(), reads.neg_reads.begin_of(chr),
                    reads.neg_reads.end_of(chr));
            setWiggleName(newfile);

            ofstream ofs(newfile.c_str());
            utils::Stamp::citationAndDate(ofs);
            rt_assert_msg(ofs.is_open(), "file not good") LOG_DEBUG1("Generating: "<<newfile);
            export_wiggle_neg(nreads, chr, ofs);
            ofs.close();
        }

    }

}

void stranded_jtwiggle::export_wiggle_gzip(Reads & reads, const char *file) {
    LOG_DEBUG1("stranded_jtwiggle::export_wiggle_gzip");
    string sfile(file);
    size_t ind = sfile.find_last_of(".wig");
    if (ind != string::npos) {
        //remove .wig
        sfile = sfile.substr(0, ind - 3);
    }

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
            preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
            npreadsend, nppreadsend, npppreadsend;
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads, mergedchrs);
    string pf(sfile + "_Pos.wig.gz");
    string nf(sfile + "_Neg.wig.gz");

    ogzstream pof(pf.c_str());
    ogzstream nof(nf.c_str());
    utils::Stamp::citationAndDate(pof);
    utils::Stamp::citationAndDate(nof);
    rt_assert_msg(pof.is_open(), "file not good")
    rt_assert_msg(nof.is_open(), "file not good")
    pof << "track type=wiggle_0 name=\"" << pf << "\" " << "visibility=dense "
            << "color=" << _colorRGB[0] << "," << _colorRGB[1] << ","
            << _colorRGB[2] << " " << "altColor=" << _colorRGB[0] << ","
            << _colorRGB[1] << "," << _colorRGB[2] << " " << "priority="
            << _priority << "\n";
    nof << "track type=wiggle_0 name=\"" << nf << "\" " << "visibility=dense "
            << "color=" << _ncolorRGB[0] << "," << _ncolorRGB[1] << ","
            << _ncolorRGB[2] << " " << "altColor=" << _ncolorRGB[0] << ","
            << _ncolorRGB[1] << "," << _ncolorRGB[2] << " " << "priority="
            << _priority << "\n";

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

        uint32_t noofbins = pchrlength ? 1 + (pchrlength / _binlength) : 0;
        uint32_t binind = 0;
        uint32_t binstart = 0, binend = 0;

        uint32_t nnoofbins = nchrlength ? 1 + (nchrlength / _binlength) : 0;
        uint32_t nbinind = 0;
        uint32_t nbinstart = 0, nbinend = 0;
        pof << "variableStep chrom=" << chr << " span=1\n";
        nof << "variableStep chrom=" << chr << " span=1\n";
        assert_gt(_binlength, 1)
        while (noofbins-- > 0) {
            LOG_DEBUG2("start Bin: "<<binind);
            binstart = _binlength * binind + 1;
            binend = _binlength * (binind + 1);
            binind++;

            if (preadsstart != preadsend) {
                ppreadsstart = lower_bound(preadsstart, preadsend, binstart);
                ppreadsend = upper_bound(ppreadsstart, preadsend, binend);
            } else {
                ppreadsstart = ppreadsend;
            }

            nppreadsstart = nppreadsend;

            JT_wiggle_file::_process(binstart, binend, reads.getReadlength(),
                    _readextlength, ppreadsstart, ppreadsend, nppreadsstart,
                    nppreadsend, pof);
            LOG_DEBUG2("pos Bin processed: "<<binind);

        }

        LOG_DEBUG1("In export_wiggle, finished processing chr pos: "<<chr);
        while (nnoofbins-- > 0) {
            LOG_DEBUG2("start Bin: "<<nbinind);
            nbinstart = _binlength * nbinind + 1;
            nbinend = _binlength * (nbinind + 1);
            nbinind++;

            if (npreadsstart != npreadsend) {
                nppreadsstart = lower_bound(npreadsstart, npreadsend,
                        nbinstart);
                nppreadsend = upper_bound(nppreadsstart, npreadsend, nbinend);
            } else {
                nppreadsstart = nppreadsend;
            }
            _process_neg(nbinstart, nbinend, reads.getReadlength(),
                    _readextlength, nppreadsstart, nppreadsend, nof);
            LOG_DEBUG2("neg Bin processed: "<<binind);

        }

        LOG_DEBUG1("In export_wiggle, finished processing chr neg: "<<chr);
    }

    LOG_DEBUG1("QUIT:stranded_jtwiggle::export_wiggle_gzip");
}

void stranded_jtwiggle::export_wiggle_pos(vector<uint32_t> & preads, string chr,
        ostream & os) {
    //    cout << "In jT_sp_export_wiggle" << endl;
    assert_neq(chr, "")
    os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
            << "visibility=dense " << "color=" << _colorRGB[0] << ","
            << _colorRGB[1] << "," << _colorRGB[2] << " " << "altColor="
            << _colorRGB[0] << "," << _colorRGB[1] << "," << _colorRGB[2] << " "
            << "priority=" << _priority << "\n";

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
            preadsend, ppreadsend, pppreadsend, nppreadsstart, nppreadsend;

    LOG_DEBUG1("In stranded_jtwiggle::export_wiggle_pos, start processing chr: "<<chr);
    preadsstart = preads.begin();
    preadsend = preads.end();

    /*
     * in case some chrs only contain pos or neg reads
     * do not use
     * pchrlength = (*(preadsend-1));
     */
    //ignore the strand if it contains less than 2 read
    //todo: this can not rule out the case preadsend = 0x00 if this function
    // is not called after reads correction.
    uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));

    uint32_t apchrlength = pchrlength;
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
        nppreadsstart = nppreadsend;
        JT_wiggle_file::_process(binstart, binend, _readlength, _readextlength,
                ppreadsstart, ppreadsend, nppreadsstart, nppreadsend, os);
        LOG_DEBUG2("Bin processed: "<<binind);

    }LOG_DEBUG1("In export_wiggle, finished processing chr: "<<chr);

}

void stranded_jtwiggle::export_wiggle_neg(vector<uint32_t> & nreads, string chr,
        ostream & os) {
    //    cout << "In jT_sp_export_wiggle" << endl;
    assert_neq(chr, "")
    os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
            << "visibility=dense " << "color=" << _ncolorRGB[0] << ","
            << _ncolorRGB[1] << "," << _ncolorRGB[2] << " " << "altColor="
            << _ncolorRGB[0] << "," << _ncolorRGB[1] << "," << _ncolorRGB[2]
            << " " << "priority=" << _priority << "\n";

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
            npreadsend, nppreadsend, npppreadsend;

    LOG_DEBUG1("In stranded_jtwiggle::export_wiggle_neg, start processing chr: "<<chr);

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
    uint32_t nchrlength = npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));

    uint32_t apchrlength = nchrlength;
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

        nppreadsstart = lower_bound(npreadsstart, npreadsend, binstart);
        nppreadsend = upper_bound(nppreadsstart, npreadsend, binend);

        this->_process_neg(binstart, binend, _readlength, _readextlength,
                nppreadsstart, nppreadsend, os);
        LOG_DEBUG2("Bin processed: "<<binind);

    }LOG_DEBUG1("In export_wiggle, finished processing chr: "<<chr);

}

void stranded_jtwiggle::split_export_wiggle_gzip(Reads & reads,
        const char *file) {
    //    cout << "In master jT_sp_export_wiggle" << endl;
    LOG_DEBUG1("stranded_jtwiggle::split_export_wiggle");
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

        if (reads.pos_reads.hasReadsOnChr(chr)) {
            vector<uint32_t> preads;

            string newfile(sfile + "_"+chr + "_Pos.wig.gz");
            preads.insert(preads.begin(), reads.pos_reads.begin_of(chr),
                    reads.pos_reads.end_of(chr));

            setWiggleName(newfile);

            ogzstream ofs(newfile.c_str());
            utils::Stamp::citationAndDate(ofs);
            rt_assert_msg(ofs.is_open(), "file not good") LOG_DEBUG1("Generating: "<<newfile);
            export_wiggle_pos(preads, chr, ofs);
            ofs.close();
        }

        if (reads.neg_reads.hasReadsOnChr(chr)) {

            vector<uint32_t> nreads;
            string newfile(sfile + "_"+chr + "_Neg.wig.gz");
            nreads.insert(nreads.begin(), reads.neg_reads.begin_of(chr),
                    reads.neg_reads.end_of(chr));
            setWiggleName(newfile);

            ogzstream ofs(newfile.c_str());
            utils::Stamp::citationAndDate(ofs);
            rt_assert_msg(ofs.is_open(), "file not good") LOG_DEBUG1("Generating: "<<newfile);
            export_wiggle_neg(nreads, chr, ofs);
            ofs.close();
        }

    }
}

//void stranded_jtwiggle::export_wiggle(Reads & reads,
//                                      data_type & os) {

//    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
//    preadsend, ppreadsend, pppreadsend;
//
//    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
//    npreadsend, nppreadsend, npppreadsend;
//    vector<string> mergedchrs;
//    reads_tools::get_merged_chrs_for_both_strands(reads,
//                                                  mergedchrs);
//    profile_zoom zm;
//    data_type tp;
//    uint32_t ws = 100000; // 100kb
//    uint32_t overlap = ws - 10000;
//    foreach(string chr, mergedchrs)
//    {
//        cout << chr << endl;
//        LOG_DEBUG1("In export_wiggle, start processing chr: "<<chr);
//        if (reads.pos_reads.hasReadsOnChr(chr)) {
//            preadsstart = reads.pos_reads.begin_of(chr);
//            preadsend = reads.pos_reads.end_of(chr);
//        } else {
//            preadsstart = preadsend;
//        }
//        if (reads.neg_reads.hasReadsOnChr(chr)) {
//            npreadsstart = reads.neg_reads.begin_of(chr);
//            npreadsend = reads.neg_reads.end_of(chr);
//        } else {
//            npreadsend = npreadsstart;
//        }
//
//        /*
//         * in case some chrs only contain pos or neg reads
//         * do not use
//         * pchrlength = (*(preadsend-1));
//         */
//        //ignore the strand if it contains less than 2 read
//        //todo: this can not rule out the case preadsend = 0x00 if this function
//        // is not called after reads correction.
//        uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
//        uint32_t nchrlength =
//        npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));
//
//        uint32_t noofbins = pchrlength ? 1 + (pchrlength / _binlength) : 0;
//        uint32_t binind = 0;
//        uint32_t binstart = 0, binend = 0;
//
//        uint32_t nnoofbins = nchrlength ? 1 + (nchrlength / _binlength) : 0;
//        uint32_t nbinind = 0;
//        uint32_t nbinstart = 0, nbinend = 0;
//
//        assert_gt(_binlength,
//                  1)
//        while (noofbins-- > 0) {
//            LOG_DEBUG2("start Bin: "<<binind);
//            binstart = _binlength * binind + 1;
//            binend = _binlength * (binind + 1);
//            binind++;
//
//            if (preadsstart != preadsend) {
//                ppreadsstart = lower_bound(preadsstart,
//                                           preadsend,
//                                           binstart);
//                ppreadsend = upper_bound(ppreadsstart,
//                                         preadsend,
//                                         binend);
//            } else {
//                ppreadsstart = ppreadsend;
//            }
//
//            nppreadsstart = nppreadsend;
//
//            _process(binstart,
//                     binend,
//                     reads.getReadlength(),
//                     _readextlength,
//                     ppreadsstart,
//                     ppreadsend,
//                     tp);
//            LOG_DEBUG2("pos Bin processed: "<<binind);
//
//        }
//        zm.smooth(ws,
//                    overlap,
//                    tp,
//                    os);
//        LOG_DEBUG1("In export_wiggle, finished processing chr pos: "<<chr);
////        while (nnoofbins-- > 0) {
////            LOG_DEBUG2("start Bin: "<<nbinind);
////            nbinstart = _binlength * nbinind + 1;
////            nbinend = _binlength * (nbinind + 1);
////            nbinind++;
////
////            ppreadsstart = ppreadsend;
////
////            if (npreadsstart != npreadsend) {
////                nppreadsstart = lower_bound(npreadsstart,
////                                            npreadsend,
////                                            nbinstart);
////                nppreadsend = upper_bound(nppreadsstart,
////                                          npreadsend,
////                                          nbinend);
////            } else {
////                nppreadsstart = nppreadsend;
////            }
////            _process_neg(nbinstart,
////                         nbinend,
////                         reads.getReadlength(),
////                         _readextlength,
////                         nppreadsstart,
////                         nppreadsend,
////                         nof);
////            LOG_DEBUG2("neg Bin processed: "<<nbinind);
////
////        }
//
//        LOG_DEBUG1("In export_wiggle, finished processing chr neg: "<<chr);
//    }

//}

void stranded_jtwiggle::_process_neg(uint32_t start, uint32_t end,
        uint32_t readlength, uint32_t readextlength,
        vector<uint32_t>::iterator nreadsStart,
        vector<uint32_t>::iterator nreadsEnd, ostream & os) {
    //    cout << "in JT_process" << endl;
    assert_gt(end, 2)
    assert_gt(end - 2, start)
    assert_gt(end - start + 1, readlength)
    assert_gt(end - start + 1, readextlength) LOG_DEBUG1("stranded_jtwiggle::_process_neg");LOG_DEBUG1("readlength:"<<readlength);LOG_DEBUG1("readextlength:"<<readextlength);
    uint32_t a;
    uint32_t b;
    uint32_t arrayStart;
    uint32_t arrayEnd;
    uint32_t read;
    typedef vector<pair<uint32_t, double> > reads_count_t;
    reads_count_t reads_count;

    bool inRange = false;

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
        //        cout <<"ps score:"<<(int32_t)(reads_count[i].second)<<"\tps pos:"<<pos<<endl;
        if (pos == ppos) {
            count += (reads_count[i].second);
            continue;
        }
        if (ppos > 0 && count > 0) {
            //This info is viewable in the result wig file
            os << ppos << "\t-" << count << "\n";
        }
        count += (reads_count[i].second);
        ppos = pos;
    }

    if (ppos > 0 && count > 0) {
        os << ppos << "\t-" << count << "\n";
    }

    LOG_DEBUG1("QUIT: stranded_jtwiggle::_process_neg");
}
//
//void stranded_jtwiggle::_process(uint32_t start,
//                                 uint32_t end,
//                                 uint32_t readlength,
//                                 uint32_t readextlength,
//                                 vector<uint32_t>::iterator nreadsStart,
//                                 vector<uint32_t>::iterator nreadsEnd,
//                                 data_type & os) {
//
////    assert_gt(end,
////              2)
////
////    assert_gt(end - 2,
////              start)
////
////    assert_gt(end - start + 1,
////              readlength)
////
////    assert_gt(end - start + 1,
////              readextlength)
////
////    LOG_DEBUG1("stranded_jtwiggle::_process_neg");
////
////    LOG_DEBUG1("readlength:"<<readlength);
////
////    LOG_DEBUG1("readextlength:"<<readextlength);
////    uint32_t a;
////    uint32_t b;
////    uint32_t arrayStart;
////    uint32_t arrayEnd;
////    uint32_t read;
////    typedef data_type reads_count_t;
////    reads_count_t reads_count;
////
////    bool inRange = false;
////
////    LOG_DEBUG1("start mapping neg reads");
////
////    LOG_DEBUG1("Total neg reads: "<<nreadsEnd-nreadsStart);
////    while (nreadsStart != nreadsEnd) {
////
////        read = *nreadsStart;
////        nreadsStart++;
////        if (read + readlength < readextlength) continue;
////        // a = read + readlength - readextlength;
////        //        b = read;
////        //THIS IS WHERE JT style kicks in
////
////        a = read + readlength - readextlength;
////        b = read + readlength;
////        //        cout <<"wig negread:"<<read<< " a:"<<a<<" b:"<<b<<endl;
////        arrayStart = 0;
////        arrayEnd = 0;
////
////        if (a < end && b > start) {
////            inRange = true;
////            /*
////             *     |-------|
////             *   |---|
////             */
////            if (a <= start && b <= end) {
////
////                reads_count.push_back(pair<uint32_t, double>(start,
////                                                              1));
////                reads_count.push_back(pair<uint32_t, double>(b,
////                                                              -1));
////            }
////            /*
////             *     |-------|
////             *       |---|
////             */
////            if (a >= start && b <= end) {
////
////                reads_count.push_back(pair<uint32_t, double>(a,
////                                                              1));
////                reads_count.push_back(pair<uint32_t, double>(b,
////                                                              -1));
////            }
////            /*
////             *     |-------|
////             *            |---|
////             */
////            if (a < end && b > end) {
////
////                reads_count.push_back(pair<uint32_t, double>(a,
////                                                              1));
////                reads_count.push_back(pair<uint32_t, double>(end,
////                                                              1));
////            }
////            /*
////             *     |-------|
////             *   |-----------|
////             */
////            if (a < start && b > end) {
////
////            }
////
////        } else if (inRange) {
////            break;
////        }
////    }
////    if (reads_count.size() == 0) {
////        return;
////    }
////    if (reads_count.size())
////    //sort based on their locations;
////    sort(reads_count.begin(),
////         reads_count.end(),
////         wiggle_sgr_sort_comparator);
////
////    uint32_t pos;
////    uint32_t ppos = (reads_count.begin())->first;
////    double count = (reads_count.begin())->second;
////    LOG_DEBUG1("start building sgr");LOG_DEBUG1("initial count value:"<<count);
////    element_type _e;
////    for (size_t i = 1; i < reads_count.size(); i++) {
////        pos = reads_count[i].first;
////
////        if (pos == ppos) {
////            count +=  (reads_count[i].second);
////            continue;
////        }
////        if (ppos > 0 && count > 0) {
////
////            _e.first = ppos;
////            _e.second = count;
////            os.push_back(_e);
////        }
////        count += (reads_count[i].second);
////        ppos = pos;
////    }
////
////    if (ppos > 0 && count > 0) {
////        _e.first = ppos;
////        _e.second = count;
////        os.push_back(_e);
////    }
////
////    LOG_DEBUG1("QUIT: stranded_jtwiggle::_process_neg");
//}
//
