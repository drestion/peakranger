/*
 * testwig.cpp
 *
 *  Created on: Jan 13, 2012
 *      Author: xfeng
 */

#include "testwig.h"
#include "wigbuilder.h"
#include "short_reads/readstools.h"
#include "wiggle/gzstream.h"

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
void test_wig::_print_wigfile_trackheader(ostream & pof,
                                          string& _name,
                                          vector<uint32_t> col)
                                          {
    pof << "track type=wiggle_0 name=\"" << _name << "\" " << "visibility=dense " << "color="
    << col[0] << "," << col[1] << "," << col[2] << " "
    << "altColor=" << col[0] << "," << col[1] << "," << col[2]
    << " " << "priority=" << _priority << "\n";
}

void test_wig::export_wiggle(Reads & reads,
                             ostream & pof,
                             string& pf) {

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
    preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
    npreadsend, nppreadsend, npppreadsend;
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads,
                                                  mergedchrs);

    rt_assert_msg(pof.good(),
                  "file not good")

    ct.print_msg(pof);

    _print_wigfile_trackheader(pof,

                               pf,
                               _colorRGB);
    wigs _wigs, _wigss;
    wig_builder _wb;
    foreach(string chr, mergedchrs)
    {
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
        pof << "variableStep chrom=" << chr << " span=1\n";
        _wb._binned_wig_compiler(_binlength,
                                 _readlength,
                                 _readextlength,
                                 preadsstart,
                                 preadsend,
                                 npreadsstart,
                                 npreadsend,
                                 pof);

    }

}

void test_wig::export_wiggle(vector<uint32_t> & preads,
                             vector<uint32_t> & nreads,
                             string chr,
                             ostream & os) {
}

//void test_wig::export_wiggle(Reads & reads,
//                             const char *file,
//                             bool gzip) {
//
//    string sfile(file);
//
//    size_t ind = sfile.find_last_of(".wig");
//    if (ind != string::npos) {
//        //remove .wig
//        sfile = sfile.substr(0,
//                             ind - 3);
//    }
//
//    string pf(sfile + ".wig");
//    ofstream pof(pf.c_str());
//    ind = sfile.find_last_of("/\\");
//    pf = sfile.substr(ind + 1);
//
//    export_wiggle(reads,
//                  pof,
//                  pf);
//
//    pof.close();
//
//}

void test_wig::_compile_neg_strand(Reads & reads,
                                   string& chr,
                                   vector<uint32_t>::iterator npreadsstart,
                                   vector<uint32_t>::iterator npreadsend,
                                   ofstream & pof)
                                   {
    if (reads.neg_reads.hasReadsOnChr(chr)) {
        npreadsstart = reads.neg_reads.begin_of(chr);
        npreadsend = reads.neg_reads.end_of(chr);
    } else {
        npreadsend = npreadsstart;
    }
    pof << "variableStep chrom=" << chr << " span=1\n";
    _w._binned_wig_compiler(_binlength,
                            _readlength,
                            _readextlength,
                            npreadsstart,
                            npreadsend,
                            pof,
                            "-",
                            wig_builder::_get_ab_re);
}

void test_wig::_compile_pos_strand(Reads & reads,
                                   string& chr,
                                   vector<uint32_t>::iterator npreadsstart,
                                   vector<uint32_t>::iterator npreadsend,
                                   ofstream & pof)
                                   {
    if (reads.pos_reads.hasReadsOnChr(chr)) {
        npreadsstart = reads.pos_reads.begin_of(chr);
        npreadsend = reads.pos_reads.end_of(chr);
    } else {
        npreadsend = npreadsstart;
    }
    pof << "variableStep chrom=" << chr << " span=1\n";
    _w._binned_wig_compiler(_binlength,
                            _readlength,
                            _readextlength,
                            npreadsstart,
                            npreadsend,
                            pof,
                            "",
                            wig_builder::_get_ab);
}
//void test_wig::export_wiggle_neg(Reads & reads,
//                                 ostream& pof,
//                                 string& pf) {
//    vector<string> mergedchrs;
//    reads_tools::get_merged_chrs_for_both_strands(reads,
//                                                  mergedchrs);
//
//    vector<uint32_t>::iterator npreadsstart, npreadsend;
//    rt_assert_msg(pof.good(),
//                  "file not good")
//
//    ct.print_msg(pof);
//
//    _print_wigfile_trackheader(pof,
//                               pf,
//                               _ncolorRGB);
//
//    foreach(string chr, mergedchrs)
//    {
//        LOG_DEBUG1("Processing "<<chr);
//        _compile_neg_strand(reads,
//                            chr,
//                            npreadsstart,
//                            npreadsend,
//                            pof
//                            );
//
//    }
//
//}
//
//void test_wig::export_wiggle_neg(Reads & reads,
//                                 const char *file,
//                                 bool gzip) {
//
//    string sfile(file);
//
//    size_t ind = sfile.find_last_of(".wig");
//    if (ind != string::npos) {
//        //remove .wig
//        sfile = sfile.substr(0,
//                             ind - 3);
//    }
//
//    string pf(sfile + "_neg.wig");
//    ind = sfile.find_last_of("/\\");
//    pf = sfile.substr(ind + 1);
//    if (gzip) {
//        ogzstream pof(pf.c_str());
//
//        export_wiggle_neg(reads,
//                          pof,
//                          pf);
//
//        pof.close();
//    }
//    else {
//        ofstream pof(pf.c_str());
//
//        export_wiggle_neg(reads,
//                          pof,
//                          pf);
//
//        pof.close();
//    }
//
//}

//export template<typename _Trans>
//void test_wig::export_wiggle_strand(vector<uint32_t>::iterator _rs,
//                                    vector<uint32_t>::iterator _re,
//                                    ostream& os,
//                                    const char* neg,
//                                    _Trans trans) {
//
//    _w._binned_wig_compiler(_binlength,
//                            _readlength,
//                            _readextlength,
//                            _rs,
//                            _re,
//                            os,
//                            neg,
//                            trans);
//
//}
void test_wig::split_export_wiggle(Reads & reads,
                                   ostream & os) {

}
//void test_wig::export_wiggle_stranded(Reads & reads,
//                                      const char *file,
//                                      bool gzip) {
//    string sfile(file);
//
//    size_t ind = sfile.find_last_of(".wig");
//    if (ind != string::npos) {
//        //remove .wig
//        sfile = sfile.substr(0,
//                             ind - 3);
//    }
//
//    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
//    preadsend, ppreadsend, pppreadsend;
//
//    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
//    npreadsend, nppreadsend, npppreadsend;
//    vector<string> mergedchrs;
//    reads_tools::get_merged_chrs_for_both_strands(reads,
//                                                  mergedchrs);
//
//    rt_assert_msg(pof.good(),
//                  "file not good")
//
//    ct.print_msg(pof);
//
//    _print_wigfile_trackheader(pof,
//                               pf);
//    wigs _wigs, _wigss;
//
//    foreach(string chr, mergedchrs)
//    {
//
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
//        if (gzip) {
//            string pf(sfile + chr + "_pos.wig");
//            ogzstream pof(pf.c_str());
//            ind = sfile.find_last_of("/\\");
//            pf = sfile.substr(ind + 1);
//            pof << "variableStep chrom=" << chr << " span=1\n";
//            _w._binned_wig_compiler(_binlength,
//                                    _readlength,
//                                    _readextlength,
//                                    npreadsstart,
//                                    npreadsend,
//                                    pof,
//                                    "-",
//                                    wig_builder::_get_ab_re
//                                    );
//        }
//        else {
//            string pf(sfile + chr + "_pos.wig");
//            ofstream pof(pf.c_str());
//            ind = sfile.find_last_of("/\\");
//            pf = sfile.substr(ind + 1);
//            pof << "variableStep chrom=" << chr << " span=1\n";
//            _w._binned_wig_compiler(_binlength,
//                                    _readlength,
//                                    _readextlength,
//                                    npreadsstart,
//                                    npreadsend,
//                                    pof,
//                                    "-",
//                                    wig_builder::_get_ab_re
//                                    );
//        }
//
//    }
//    pof.close();
//}
//void test_wig::split_export_wiggle(Reads & reads,
//                                   const char *file,
//                                   bool gzip,
//                                   bool stranded) {
//    //    cout << "In master jT_sp_export_wiggle" << endl;
//    LOG_DEBUG1("JT_wiggle_file::split_export_wiggle");
//    rt_assert(file)
//    string sfile(file);
//    size_t ind = sfile.find_last_of(".wig");
//    if (ind != string::npos) {
//        //remove .wig
//        sfile = sfile.substr(0,
//                             ind - 3);
//    }
//    vector<string> mergedchrs;
//    reads_tools::get_merged_chrs_for_both_strands(reads,
//                                                  mergedchrs);
//
//    foreach(string chr, mergedchrs)
//    {
//        setWiggleName(sfile + chr + ".wig");
//
//        string newfile(sfile + chr + ".wig");
//        ofstream ofs(newfile.c_str());
//        rt_assert_msg(ofs.is_open(),
//                      "file not good")
//        vector<uint32_t> preads;
//
//        vector<uint32_t> nreads;
//        if (reads.pos_reads.hasReadsOnChr(chr)) {
//            preads.insert(preads.begin(),
//                          reads.pos_reads.begin_of(chr),
//                          reads.pos_reads.end_of(chr));
//        }
//
//        if (reads.neg_reads.hasReadsOnChr(chr)) {
//            nreads.insert(nreads.begin(),
//                          reads.neg_reads.begin_of(chr),
//                          reads.neg_reads.end_of(chr));
//        }
//
//        LOG_DEBUG1("Generating: "<<newfile);
//        export_wiggle(preads,
//                      nreads,
//                      chr,
//                      ofs);
//        ofs.close();
//    }
//}

void test_wig::export_wiggle_gzip(Reads & reads,
                                  const char *file) {
}

void test_wig::split_export_wiggle_gzip(Reads & reads,
                                        const char *file) {
}

