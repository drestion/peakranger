/*
 * calledpeak.h
 *
 *  Created on: Dec 9, 2011
 *      Author: xfeng
 */

#ifndef CALLEDPEAK_H_
#define CALLEDPEAK_H_

#include <stdint.h>
#include <vector>
#include <iostream>

class called_peak {
    friend std::ostream& operator<<(std::ostream& os, const called_peak& pk);
public:

    called_peak() {
        first = 0;
        second = 0;
        p = 0;
        q = 0;
        summits.resize(0);
        treads = 0;
        creads = 0;
    }

    called_peak(uint32_t f, uint32_t s, double _p, double fd, uint32_t tr,
            uint32_t cr, std::vector<uint32_t> sum) {
        first = f;
        second = s;
        p = _p;
        q = fd;
        treads = tr;
        creads = cr;
        summits.resize(0);
        summits.insert(summits.begin(), sum.begin(), sum.end());
    }

    uint32_t first;
    uint32_t second;
    double p;
    double q;
    uint32_t treads;
    uint32_t creads;
    std::vector<uint32_t> summits;

};

#endif /* CALLEDPEAK_H_ */
