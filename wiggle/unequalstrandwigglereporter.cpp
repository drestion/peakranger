/*
 * unequalstrandwigglereporter.cpp
 *
 *  Created on: Sep 27, 2011
 *      Author: xfeng
 */

#include "unequalstrandwigglereporter.h"
#include "utils/logger.h"
#include "utils/exceptions.h"
#include "utils/assert_helpers.h"
#include "short_reads/readstools.h"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include <math.h>
#include <stdint.h>
#include <algorithm>
#include <fstream>
#include <memory>
#include <vector>
#include <sstream>

using namespace std;
#define foreach BOOST_FOREACH

unequal_strand_wiggle_reporter::unequal_strand_wiggle_reporter()
{
    
}

unequal_strand_wiggle_reporter::~unequal_strand_wiggle_reporter()
{

}

void unequal_strand_wiggle_reporter::export_wiggle(Reads & reads,
                                                   ostream & os)
                                                   {

    os << "track type=wiggle_0 name=\"" << getWiggleName() << "\" "
    << "visibility=dense autoScale=off " << "viewLimits=" << _viewLimitDown
    << ":" << _viewLimitUp << " " << "color=" << _colorRGB[0] << ","
    << _colorRGB[1] << "," << _colorRGB[2] << " " << "altColor="
    << _colorRGB[0] << "," << _colorRGB[1] << "," << _colorRGB[2] << " "
    << "priority=" << _priority << "\n";

    vector<uint32_t>::iterator preadsstart, ppreadsstart, pppreadsstart,
    preadsend, ppreadsend, pppreadsend;

    vector<uint32_t>::iterator npreadsstart, nppreadsstart, npppreadsstart,
    npreadsend, nppreadsend, npppreadsend;
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads,
                                                  mergedchrs);

    foreach(string chr, mergedchrs)
    {
        LOG_DEBUG1("In export_wiggle, start processing chr: "<<chr);
        if (reads.pos_reads.hasReadsOnChr(chr)) {
            preadsstart = reads.pos_reads.begin_of(chr);
            preadsend = reads.pos_reads.end_of(chr);
        }
        else {
            preadsstart = preadsend;
        }
        if (reads.neg_reads.hasReadsOnChr(chr)) {
            npreadsstart = reads.neg_reads.begin_of(chr);
            npreadsend = reads.neg_reads.end_of(chr);
        }
        else {
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
        uint32_t nchrlength = npreadsend == npreadsstart ? 0 : (*(npreadsend
        - 1));

        uint32_t apchrlength = pchrlength > nchrlength ? pchrlength
        : nchrlength;
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

            if (preadsstart != preadsend) {
                ppreadsstart = lower_bound(preadsstart,
                                           preadsend,
                                           binstart);
                ppreadsend = upper_bound(ppreadsstart,
                                         preadsend,
                                         binend);
            }
            else {
                ppreadsstart = ppreadsend;
            }
            if (npreadsstart != npreadsend) {
                nppreadsstart = lower_bound(npreadsstart,
                                            npreadsend,
                                            binstart);
                nppreadsend = upper_bound(nppreadsstart,
                                          npreadsend,
                                          binend);
            }
            else {
                nppreadsstart = nppreadsend;
            }
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

void unequal_strand_wiggle_reporter::split_export_wiggle(Reads & reads,
                                                         const char *file)
                                                         {
    LOG_DEBUG1("unequal_strand_wiggle_reporter::split_export_wiggle");
    rt_assert(file)
    string sfile(file);
    size_t ind = sfile.find_last_of(".wig");
    if (ind != string::npos) {
        //remove .wig
        sfile = sfile.substr(0,
                             ind - 3);
    }
    vector<string> mergedchrs;
    reads_tools::get_merged_chrs_for_both_strands(reads,
                                                  mergedchrs);
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
    foreach(string chr, mergedchrs)
    {
        setWiggleName(sfile + chr + ".wig");

        string newfile(sfile + chr + ".wig");
        ofstream ofs(newfile.c_str());
        rt_assert_msg(ofs.is_open(),"file not good")
        vector<uint32_t> preads;

        vector<uint32_t> nreads;
        if (reads.pos_reads.hasReadsOnChr(chr)) {
            preads.insert(preads.begin(),
                          reads.pos_reads.begin_of(chr),
                          reads.pos_reads.end_of(chr));
        }

        if (reads.neg_reads.hasReadsOnChr(chr)) {
            nreads.insert(nreads.begin(),
                          reads.neg_reads.begin_of(chr),
                          reads.neg_reads.end_of(chr));
        }

        LOG_DEBUG1("Generating: "<<newfile);
        wiggle_reporter::export_wiggle(preads,
                                       nreads,
                                       chr,
                                       ofs);
        ofs.close();
    }
}

void unequal_strand_wiggle_reporter::export_wiggle(Reads & reads,
                                                   const char *file)
                                                   {

    setWiggleName(string(file));
    use_default_setting();
    ofstream ofs(file);
    rt_assert_msg(ofs.is_open(),"file not good")
    export_wiggle(reads,
                  ofs);
    ofs.close();
}

