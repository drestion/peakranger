/*
 * chip4cwig.h
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#ifndef CHIP4CWIG_H_
#define CHIP4CWIG_H_

#include "zoomable.h"
#include "short_reads/reads.h"
#include "wigbuilder.h"
#include "utils/util_print.h"
#include "headerprintable.h"
#include "containsdefaultposwigcolor.h"

#include <stdint.h>

class chip_4c_wig: public zoomable, public header_printable, public contains_default_pos_wig_color {
public:
    chip_4c_wig();
    virtual ~chip_4c_wig();

    virtual void smooth(wigs& _wigs, wigs& result, uint32_t window, uint32_t overlap);
    void export_smoothed_wig(Reads& reads, const char * file, bool gzip = true);
    void export_splitted_smoothed_wig(Reads& reads, const char * file, bool gzip = true);
    uint32_t getBinlength() const;
    uint32_t getOv() const;
    uint32_t getReadextlength() const;
    uint32_t getReadlength() const;
    wig_builder getW() const;
    uint32_t getWs() const;
    void setBinlength(uint32_t _binlength);
    void setOv(uint32_t _ov);
    void setReadextlength(uint32_t _readextlength);
    void setReadlength(uint32_t _readlength);
    void setW(wig_builder _w);
    void setWs(uint32_t _ws);
    uint32_t getRepeat() const;
    void setRepeat(uint32_t _repeat);

protected:
    wig_builder _w;
    utilprint::citation _ct;
    uint32_t _ws;
    uint32_t _ov;
    uint32_t _readextlength;
    uint32_t _readlength;
    uint32_t _binlength;
    uint32_t _repeat;
};

#endif /* CHIP4CWIG_H_ */
