/*
 * wigbuilder.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: xfeng
 */

#include "wigbuilder.h"

#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <ostream>
using namespace std;

void wig_builder::_process(uint32_t start,
                           uint32_t end,
                           uint32_t readlength,
                           uint32_t readextlength,
                           std::vector<uint32_t>::iterator readsStart,
                           std::vector<uint32_t>::iterator readsEnd,
                           std::vector<uint32_t>::iterator nreadsStart,
                           std::vector<uint32_t>::iterator nreadsEnd,
                           wigs& r) {
    assert_gt(end,
              2)

    assert_gt(end - 2,
              start)

    assert_gt(end-start+1,
              readlength)

    assert_gt(end-start+1,
              readextlength)

    LOG_DEBUG1("Entering wig_builder::_process");

    LOG_DEBUG1("readlength:"<<readlength);

    LOG_DEBUG1("readextlength:"<<readextlength);

    _process(start,
             end,
             readlength,
             readextlength,
             readsStart,
             readsEnd,
             r,
             _get_ab);
    _process(start,
             end,
             readlength,
             readextlength,
             nreadsStart,
             nreadsEnd,
             r,
             _get_ab_re);
    sort(r.begin(),
         r.end(),
         wig::compa);
    LOG_DEBUG1("QUIT: wig_builder::_process");
}

void wig_builder::_process_wig(uint32_t start,
                               uint32_t end,
                               uint32_t readlength,
                               uint32_t readextlength,
                               std::vector<uint32_t>::iterator readsStart,
                               std::vector<uint32_t>::iterator readsEnd,
                               std::vector<uint32_t>::iterator nreadsStart,
                               std::vector<uint32_t>::iterator nreadsEnd,
                               wigs& r) {
    LOG_DEBUG1("Entering wig_builder::_process_wig");
    wigs _r;
    _process(start,
             end,
             readlength,
             readextlength,
             readsStart,
             readsEnd,
             nreadsStart,
             nreadsEnd,
             _r);
    _compile_wig(_r,
                 r);
    LOG_DEBUG1("QUIT: wig_builder::_process_wig");
}

void wig_builder::_compile_wig(wigs& _r,
                               wigs& r) {
    LOG_DEBUG1("Entering wig_builder::compile_wig");
    if (_r.size() == 0) {
        return;
    }
    //sort based on their locations;
    sort(_r.begin(),
         _r.end(),
         wig::compa);

    uint32_t pos;
    uint32_t ppos = (_r.begin())->getP();
    double count = (_r.begin())->getS();

    LOG_DEBUG1("initial count value:"<<count);
    for (size_t i = 1; i < _r.size(); i++) {
        LOG_DEBUG2("get pos:"<<pos<<" count:"<<count);
        pos = _r[i].getP();

        if (pos == ppos) {
            count += (_r[i].getS());
            continue;
        }
        if (ppos > 0 && count > 0) {
            r.push_back(wig(ppos,
                            count));
        }
        count += (_r[i].getS());
        ppos = pos;
    }
    if (ppos > 0 && count > 0) {
        r.push_back(wig(ppos,
                        count));
    }

    LOG_DEBUG1("QUIT: wig_builder::compile_wig");
}

void wig_builder::_get_ab(uint32_t read,
                          uint32_t readlength,
                          uint32_t ext,
                          uint32_t & a,
                          uint32_t & b) {
    a = read;
    b = read + ext;
}

void wig_builder::_process(uint32_t readlength,
                           uint32_t readextlength,
                           std::vector<uint32_t>::iterator readsStart,
                           std::vector<uint32_t>::iterator readsEnd,
                           std::vector<uint32_t>::iterator nreadsStart,
                           std::vector<uint32_t>::iterator nreadsEnd,
                           wigs & r) {
    LOG_DEBUG1("Entering wig_builder::_process");

    LOG_DEBUG1("readlength:"<<readlength);

    LOG_DEBUG1("readextlength:"<<readextlength);

    _process(readlength,
             readextlength,
             readsStart,
             readsEnd,
             r,
             _get_ab);
    _process(readlength,
             readextlength,
             nreadsStart,
             nreadsEnd,
             r,
             _get_ab_re);

    sort(r.begin(),
         r.end(),
         wig::compa);
}

void wig_builder::_process_wig(uint32_t readlength,
                               uint32_t readextlength,
                               std::vector<uint32_t>::iterator readsStart,
                               std::vector<uint32_t>::iterator readsEnd,
                               std::vector<uint32_t>::iterator nreadsStart,
                               std::vector<uint32_t>::iterator nreadsEnd,
                               wigs & r) {
    LOG_DEBUG1("Entering wig_builder::_process_wig");
    wigs _r;
    _process(readlength,
             readextlength,
             readsStart,
             readsEnd,
             nreadsStart,
             nreadsEnd,
             _r);
    _compile_wig(_r,
                 r);
    LOG_DEBUG1("QUIT: wig_builder::_process_wig");
}

void wig_builder::_process_wig(uint32_t start,
                               uint32_t end,
                               uint32_t readlength,
                               uint32_t readextlength,
                               std::vector<uint32_t>::iterator readsStart,
                               std::vector<uint32_t>::iterator readsEnd,
                               std::vector<uint32_t>::iterator nreadsStart,
                               std::vector<uint32_t>::iterator nreadsEnd,
                               ostream & r)
                               {
    LOG_DEBUG1("Entering wig_builder::_process_wig");
    wigs _r;
    _process(start,
             end,
             readlength,
             readextlength,
             readsStart,
             readsEnd,
             nreadsStart,
             nreadsEnd,
             _r);
    _compile_wig(_r,
                 r);
    LOG_DEBUG1("QUIT: wig_builder::_process_wig");
}
void wig_builder::_binned_wig_compiler(uint32_t _binlength,
                                       uint32_t _readlength,
                                       uint32_t _readextlength,
                                       std::vector<uint32_t>::iterator preadsstart,
                                       std::vector<uint32_t>::iterator preadsend,
                                       std::vector<uint32_t>::iterator npreadsstart,
                                       std::vector<uint32_t>::iterator npreadsend,
                                       ostream & pof)
                                       {
    LOG_DEBUG1("wig_builder::_binned_wig_compiler");
    std::vector<uint32_t>::iterator ppreadsstart, pppreadsstart,
    ppreadsend, pppreadsend, nppreadsstart, nppreadsend;
    //todo: this can not rule out the case preadsend = 0x00 if this function
    // is not called after reads correction.
    uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
    uint32_t nchrlength = npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));

    uint32_t noofbins = pchrlength ? 1 + (pchrlength / _binlength) : 0;
    uint32_t binind = 0;
    uint32_t binstart = 0, binend = 0;
    uint32_t nnoofbins = nchrlength ? 1 + (nchrlength / _binlength) : 0;
    uint32_t nbinind = 0;
    uint32_t nbinstart = 0, nbinend = 0;
    LOG_DEBUG1("pos chr length:"<<pchrlength);
    LOG_DEBUG1("neg chr length:"<<nchrlength);
    LOG_DEBUG1("pos chr bins:"<<noofbins);
    LOG_DEBUG1("pos chr bins:"<<nnoofbins);
    assert_gt(_binlength,
              1)
    bool _pb = false;
    bool _nb = false;
    wigs _wigs;
    while (true) {
        if (noofbins > 0) {
            binstart = _binlength * binind + 1;
            binend = _binlength * (binind + 1);
            LOG_DEBUG2("in pos chr bin:"<<binstart<<":"<<binend);
            binind++;
            ppreadsstart = lower_bound(preadsstart,
                                       preadsend,
                                       binstart);
            ppreadsend = upper_bound(ppreadsstart,
                                     preadsend,
                                     binend);
            _process(binstart,
                     binend,
                     _readlength,
                     _readextlength,
                     ppreadsstart,
                     ppreadsend,
                     _wigs,
                     _get_ab);

            noofbins--;
        } else {
            LOG_DEBUG2("finished all pos bins");
            _pb = true;
        }
        if (nnoofbins > 0) {
            nbinstart = _binlength * nbinind + 1;
            nbinend = _binlength * (nbinind + 1);
            LOG_DEBUG2("in neg chr bin:"<<nbinstart<<":"<<nbinend);
            nbinind++;
            nppreadsstart = lower_bound(npreadsstart,
                                        npreadsend,
                                        nbinstart);
            nppreadsend = upper_bound(nppreadsstart,
                                      npreadsend,
                                      nbinend);
            _process(nbinstart,
                     nbinend,
                     _readlength,
                     _readextlength,
                     nppreadsstart,
                     nppreadsend,
                     _wigs,
                     _get_ab_re);

            nnoofbins--;
        } else {
            LOG_DEBUG2("finished all neg bins");
            _nb = true;
        }
        _compile_wig(_wigs,
                     pof);
        _wigs.clear();
        if (_pb && _nb) {
            break;
        }
    }
}






void wig_builder::_binned_wig_compiler(uint32_t _binlength,
                                       uint32_t _readlength,
                                       uint32_t _readextlength,
                                       std::vector<uint32_t>::iterator preadsstart,
                                       std::vector<uint32_t>::iterator preadsend,
                                       std::vector<uint32_t>::iterator npreadsstart,
                                       std::vector<uint32_t>::iterator npreadsend,
                                       wigs & pof)
                                       {
    LOG_DEBUG1("wig_builder::_binned_wig_compiler");
    std::vector<uint32_t>::iterator ppreadsstart, pppreadsstart,
    ppreadsend, pppreadsend, nppreadsstart, nppreadsend;
    //todo: this can not rule out the case preadsend = 0x00 if this function
    // is not called after reads correction.
    uint32_t pchrlength = preadsend == preadsstart ? 0 : (*(preadsend - 1));
    uint32_t nchrlength = npreadsend == npreadsstart ? 0 : (*(npreadsend - 1));

    uint32_t noofbins = pchrlength ? 1 + (pchrlength / _binlength) : 0;
    uint32_t binind = 0;
    uint32_t binstart = 0, binend = 0;
    uint32_t nnoofbins = nchrlength ? 1 + (nchrlength / _binlength) : 0;
    uint32_t nbinind = 0;
    uint32_t nbinstart = 0, nbinend = 0;
    LOG_DEBUG1("pos chr length:"<<pchrlength);
    LOG_DEBUG1("neg chr length:"<<nchrlength);
    LOG_DEBUG1("pos chr bins:"<<noofbins);
    LOG_DEBUG1("pos chr bins:"<<nnoofbins);
    assert_gt(_binlength,
              1)
    bool _pb = false;
    bool _nb = false;
    wigs _wigs;
    while (true) {
        if (noofbins > 0) {
            binstart = _binlength * binind + 1;
            binend = _binlength * (binind + 1);
            LOG_DEBUG2("in pos chr bin:"<<binstart<<":"<<binend);
            binind++;
            ppreadsstart = lower_bound(preadsstart,
                                       preadsend,
                                       binstart);
            ppreadsend = upper_bound(ppreadsstart,
                                     preadsend,
                                     binend);
            _process(binstart,
                     binend,
                     _readlength,
                     _readextlength,
                     ppreadsstart,
                     ppreadsend,
                     _wigs,
                     _get_ab);

            noofbins--;
        } else {
            LOG_DEBUG2("finished all pos bins");
            _pb = true;
        }
        if (nnoofbins > 0) {
            nbinstart = _binlength * nbinind + 1;
            nbinend = _binlength * (nbinind + 1);
            LOG_DEBUG2("in neg chr bin:"<<nbinstart<<":"<<nbinend);
            nbinind++;
            nppreadsstart = lower_bound(npreadsstart,
                                        npreadsend,
                                        nbinstart);
            nppreadsend = upper_bound(nppreadsstart,
                                      npreadsend,
                                      nbinend);
            _process(nbinstart,
                     nbinend,
                     _readlength,
                     _readextlength,
                     nppreadsstart,
                     nppreadsend,
                     _wigs,
                     _get_ab_re);

            nnoofbins--;
        } else {
            LOG_DEBUG2("finished all neg bins");
            _nb = true;
        }
        _compile_wig(_wigs,
                     pof);
        _wigs.clear();
        if (_pb && _nb) {
            break;
        }
    }
}

void wig_builder::_compile_wig(wigs & _r,
                               ostream & r,
                               const char* _neg)
                               {
    LOG_DEBUG1("Entering wig_builder::compile_wig");
    if (_r.size() == 0) {
        return;
    }
    //sort based on their locations;
    sort(_r.begin(),
         _r.end(),
         wig::compa);

    uint32_t pos;
    uint32_t ppos = (_r.begin())->getP();
    double count = (_r.begin())->getS();

    LOG_DEBUG1("initial count value:"<<count);
    for (size_t i = 1; i < _r.size(); i++) {
        LOG_DEBUG2("get pos:"<<pos<<" count:"<<count);
        pos = _r[i].getP();

        if (pos == ppos) {
            count += (_r[i].getS());
            continue;
        }
        if (ppos > 0 && count > 0) {
            r << ppos << "\t" << _neg << count << "\n";
        }
        count += (_r[i].getS());
        ppos = pos;
    }
    if (ppos > 0 && count > 0) {
        r << ppos << "\t" << _neg << count << "\n";
    }

    LOG_DEBUG1("QUIT: wig_builder::compile_wig");
}


void wig_builder::_get_ab_re(uint32_t read,
                             uint32_t readlength,
                             uint32_t ext,
                             uint32_t & a,
                             uint32_t & b) {
    if (read + readlength >= ext) {
        a = read + readlength - ext;
    } else {
        a = 0;
    }
    b = read + readlength;
}

