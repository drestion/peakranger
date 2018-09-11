/*
 * point_reads.h
 *
 *  Created on: Jul 26, 2011
 *      Author: xin
 */

#ifndef POINT_READS_H_
#define POINT_READS_H_
#include "reads.h"
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>

/*
 * Mimic the peakseq style profile
 * The result is a vector of pair of pos
 */

typedef std::pair<uint32_t, int8_t> pos_t;

class point_reads {
public:
    point_reads();
    virtual ~point_reads();

    void set_reads(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd,
            std::vector<pos_t>& result);

    void set_reads(uint32_t readlength, uint32_t readextlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<pos_t>& result);

    void set_reads(std::vector<uint32_t>& reads, uint32_t extension,
            std::vector<pos_t>& result);

};

#endif /* POINT_READS_H_ */
