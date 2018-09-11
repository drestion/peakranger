/*
 * peakseq_profile.h
 *
 *  Created on: Jul 26, 2011
 *      Author: xin
 */

#ifndef PEAKSEQ_PROFILE_H_
#define PEAKSEQ_PROFILE_H_
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>

/*
 * Mimic the peakseq style profile
 * The result is a vector of pair of pos
 */

typedef std::pair<uint32_t, uint32_t> profile_pos_t;

class peakseq_profile {
public:
    peakseq_profile();
    virtual ~peakseq_profile();

    static void get_profile_of_reads(std::vector<uint32_t>& reads,
            uint32_t extension, std::vector<profile_pos_t>& result);

    static void get_profile_of_reads(uint32_t start, uint32_t end,
            uint32_t readlength, uint32_t readextlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd,
            std::vector<profile_pos_t>& result);
    static void dump(std::vector<profile_pos_t>& profile, const char* filename);

};

#endif /* PEAKSEQ_PROFILE_H_ */
