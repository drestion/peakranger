/*
 * wigbuilder.h
 *
 *  Created on: Jan 12, 2012
 *      Author: xfeng
 */

#ifndef WIGBUILDER_H_
#define WIGBUILDER_H_


#include "utils/logger.h"
#include "wig.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include "utils/stringutil.h"
#include <stdint.h>
#include <iostream>
#include <vector>

class wig_builder {
public:
    void _process_wig(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, std::ostream& r);
    void _process_wig(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, wigs& r);

    void _process_wig(uint32_t readlength, uint32_t readextlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, wigs& r);
    void _process(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, wigs& r);
    void _process(uint32_t readlength, uint32_t readextlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, wigs& r);
    template<class _T>
    void _process(uint32_t readlength, uint32_t readextlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd, wigs& r, _T trans) {

        LOG_DEBUG1("Entering wig_builder::_process");

        LOG_DEBUG1("readlength:"<<readlength);

        LOG_DEBUG1("readextlength:"<<readextlength);

        uint32_t a;
        uint32_t b;

        uint32_t read;
        while (readsStart != readsEnd) {
            read = *readsStart++;
            trans(read, readlength, readextlength, a, b);
            r.push_back(wig(a, 1));
            r.push_back(wig(b, -1));
        }

        //sort based on their locations;
        sort(r.begin(), r.end(), wig::compa);

        LOG_DEBUG1("QUIT: wig_builder::_process");
    }

    template<class _T>
    void _process(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd, wigs& r, _T trans) {

        LOG_DEBUG1("Entering wig_builder::_process");

        LOG_DEBUG1("readlength:"<<readlength);

        LOG_DEBUG1("readextlength:"<<readextlength);

        assert_gt(end, 2)

        assert_gt(end - 2, start)

        assert_gt(end-start+1, readlength)

        assert_gt(end-start+1, readextlength)

        uint32_t a;
        uint32_t b;
        uint32_t arrayStart;
        uint32_t arrayEnd;
        uint32_t read;

        bool inRange = false;
        LOG_DEBUG1("start mapping  reads");

        LOG_DEBUG1("Total reads: "<<readsEnd-readsStart);
        while (readsStart != readsEnd) {

            read = *readsStart;
            readsStart++;
            trans(read, readlength, readextlength, a, b);
            arrayStart = 0;
            arrayEnd = 0;
            LOG_DEBUG2("Get read : " <<read<<" a:"<<a <<" b:"<<b);
            if (a < end && b > start) {
                inRange = true;
                /*
                 *     |-------|
                 *   |---|
                 */
                if (a <= start && b <= end) {
                    LOG_DEBUG2("left ");
                    r.push_back(wig(start, 1));
                    r.push_back(wig(b, -1));
                }
                /*
                 *     |-------|
                 *       |---|
                 */
                if (a >= start && b <= end) {
                    LOG_DEBUG2("center ");
                    r.push_back(wig(a, 1));
                    r.push_back(wig(b, -1));
                }
                /*
                 *     |-------|
                 *            |---|
                 */
                if (a < end && b > end) {
                    LOG_DEBUG2("right ");
                    r.push_back(wig(a, 1));
                    r.push_back(wig(end, 1));
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

        if (r.size() == 0) {
            return;
        }
        //sort based on their locations;
        sort(r.begin(), r.end(), wig::compa);

        LOG_DEBUG1("QUIT: wig_builder::_process");
    }

    void _binned_wig_compiler(uint32_t binlength, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator preadsstart,
            std::vector<uint32_t>::iterator preadsend,
            std::vector<uint32_t>::iterator npreadsstart,
            std::vector<uint32_t>::iterator npreadsend, std::ostream& pof);

    void _binned_wig_compiler(uint32_t binlength, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator preadsstart,
            std::vector<uint32_t>::iterator preadsend,
            std::vector<uint32_t>::iterator npreadsstart,
            std::vector<uint32_t>::iterator npreadsend, wigs & pof);

    template<typename _T>
    void _binned_wig_compiler(uint32_t _binlength, uint32_t _readlength,
            uint32_t _readextlength, std::vector<uint32_t>::iterator preadsstart,
            std::vector<uint32_t>::iterator preadsend, std::ostream& pof,
            const char* _neg, _T trans) {
        LOG_DEBUG1("wig_builder::_binned_wig_compiler");
        std::vector<uint32_t>::iterator ppreadsstart, pppreadsstart, ppreadsend,
                pppreadsend;
        //todo: this can not rule out the case preadsend = 0x00 if this function
        // is not called after reads correction.
        uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
        uint32_t noofbins = pchrlength ? 1 + (pchrlength / _binlength) : 0;
        uint32_t binind = 0;
        uint32_t binstart = 0, binend = 0;

        LOG_DEBUG1(" chr length:"<<pchrlength);

        LOG_DEBUG1(" chr bins:"<<noofbins);

        assert_gt(_binlength, 1)
        bool _pb = false;

        wigs _wigs;
        while (true) {
            if (noofbins > 0) {
                binstart = _binlength * binind + 1;
                binend = _binlength * (binind + 1);
                LOG_DEBUG2("in bin:"<<binstart<<":"<<binend);
                binind++;
                ppreadsstart = lower_bound(preadsstart, preadsend, binstart);
                ppreadsend = upper_bound(ppreadsstart, preadsend, binend);
                _process(binstart, binend, _readlength, _readextlength,
                        ppreadsstart, ppreadsend, _wigs, trans);

                noofbins--;
            } else {
                LOG_DEBUG2("finished all bins");
                _pb = true;
            }

            _compile_wig(_wigs, pof, _neg);
            _wigs.clear();
            if (_pb) {
                break;
            }
        }
    }

    void _compile_wig(wigs& _r, wigs& r);
    void _compile_wig(wigs& _r, std::ostream& r, const char *neg = "");
    static void _get_ab(uint32_t read, uint32_t readlength, uint32_t ext,
            uint32_t& a, uint32_t& b);

    static void _get_ab_re(uint32_t read, uint32_t readlength, uint32_t ext,
            uint32_t& a, uint32_t& b);
};

#endif /* WIGBUILDER_H_ */
