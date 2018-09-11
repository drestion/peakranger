/*
 * testwig.h
 *
 *  Created on: Jan 13, 2012
 *      Author: xfeng
 */

#ifndef TESTWIG_H_
#define TESTWIG_H_
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include "short_reads/reads.h"
#include "utils/logger.h"
#include "wiggle/wigbuilder.h"
#include "utils/util_print.h"

#include <ostream>
#include <stdint.h>
#include <iostream>

class test_wig {
public:
    test_wig() {
        uint32_t ga[3] = { 50, 126, 184 };
        _colorRGB = std::vector<uint32_t>(ga, ga + 3);
        uint32_t ga2[3] = { 228, 26, 28 };
        _ncolorRGB = std::vector<uint32_t>(ga2, ga2 + 3);
        _priority = 20;
        _viewLimitDown = 0;
        _viewLimitUp = 1024;
        _binlength = 100000000;
        _readextlength = 150;
        _readlength = 25;
    }

    virtual ~test_wig() {
    }
    virtual void export_wiggle(Reads& reads, std::ostream& os) {
    }
    virtual void export_wiggle(Reads& reads, std::ostream& os, std::string& wigname);
    virtual void export_wiggle(std::vector<uint32_t>& preads,
            std::vector<uint32_t>& nreads, std::string chr, std::ostream& os);
    typedef std::vector<std::pair<uint32_t, double> > data_type;
    typedef std::pair<uint32_t, double> element_type;
    virtual void export_wiggle(Reads& reads, data_type& result) {
    }

//    virtual void export_wiggle(Reads& reads,
//                               const char* file,
//                               bool gzip=true);
    virtual void export_wiggle(Reads& reads, const char* file) {
    }
    virtual void split_export_wiggle(Reads& reads, std::ostream& os);
    virtual void split_export_wiggle(Reads& reads, const char* file) {
    }
    virtual void export_wiggle_gzip(Reads& reads, const char* file);

    virtual void split_export_wiggle_gzip(Reads& reads, const char* file);

    uint32_t _viewLimitDown;
    uint32_t _viewLimitUp;
    std::string _wiggleName;
    std::vector<uint32_t> _colorRGB;
    std::vector<uint32_t> _ncolorRGB;
    uint32_t _priority;
    uint32_t _readextlength;
    uint32_t _readlength;
    uint32_t _binlength;

    void _print_wigfile_trackheader(std::ostream & pof, std::string& pf,
            std::vector<uint32_t> col);

//    template<typename _Trans>
//    void export_wiggle_strand(std::vector<uint32_t>::iterator read_start,
//                             std::vector<uint32_t>::iterator read_end,
//                              std::ostream& os,
//                              const char* neg,
//                              _Trans neg_read_trans = wig_builder::_get_ab);

    void _compile_neg_strand(Reads & reads, std::string& chr,
            std::vector<uint32_t>::iterator npreadsstart,
            std::vector<uint32_t>::iterator npreadsend, std::ofstream & pof);

    void _compile_pos_strand(Reads & reads, std::string& chr,
            std::vector<uint32_t>::iterator npreadsstart,
            std::vector<uint32_t>::iterator npreadsend, std::ofstream & pof);
protected:
    wig_builder _w;
    utilprint::citation ct;
};

#endif /* TESTWIG_H_ */
