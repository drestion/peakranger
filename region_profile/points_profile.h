/*
 * points_profile.h
 *
 *  Created on: Jul 26, 2011
 *      Author: xin
 */

#ifndef POINTS_PROFILE_H_
#define POINTS_PROFILE_H_
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>

/*
 * This profiler should be much much
 * faster than the region_profile class.
 * It uses a different algorithm
 */
class points_profile {
public:
    points_profile();
    virtual ~points_profile();

    static void get_profile_of_reads(std::vector<uint32_t>& reads,
            uint32_t extension, std::vector<uint16_t>& result);

    static void get_profile_of_reads(uint32_t start, uint32_t end,
            uint32_t readlength, uint32_t readextlength,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd,
            std::vector<uint16_t>& result);

};

#endif /* POINTS_PROFILE_H_ */
