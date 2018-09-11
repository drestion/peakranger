/*
 * peakseq_profile_thresholder.h
 *
 *  Created on: July 26, 2011
 *      Author: xin
 */

#ifndef PEAKSEQ_PROFILE_THRESHOLDER_H_H
#define PEAKSEQ_PROFILE_THRESHOLDER_H_H

#include <vector>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <stdint.h>
#include "utils/exceptions.h"

typedef std::pair<uint32_t, uint32_t> profile_pos_t;
typedef std::pair<uint32_t, uint32_t> peak_pos_t;
class peakseq_profile_thresholder {
public:

    peakseq_profile_thresholder() {
    }
    ~peakseq_profile_thresholder() {
    }

public:
    /*
     * result will be reset
     */
    static void threshold(std::vector<profile_pos_t>& vec, uint32_t cutoff,
            uint32_t mergeDistance, uint32_t offset,
            std::vector<peak_pos_t>& result);
};

#endif /* PEAKSEQ_PROFILE_THRESHOLDER_H_H */
