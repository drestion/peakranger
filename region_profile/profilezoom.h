/*
 * profilezoom.h
 *  Designed for peakseq style profile.
 *  Created on: Jan 11, 2012
 *      Author: xfeng
 */

#ifndef PROFILEZOOM_H_
#define PROFILEZOOM_H_

#include <stdint.h>
#include <vector>
#include <utility>
#include <iostream>
#include <queue>

#include "utils/assert_helpers.h"
#include "utils/logger.h"
#include "wiggle/wig.h"
class profile_zoom {
public:

    typedef std::vector<wig> data_type;
    typedef wig element_type;

    const void zoom_out(const uint32_t window_size,
            const data_type data_to_zoom, data_type& result);

    const void zoom_out(const uint32_t window_size, const uint32_t overlap_size,
            const data_type data_to_zoom, data_type& result);

    const void smooth(const uint32_t window_size, const uint32_t overlap_size,
            const data_type data_to_zoom, data_type& result);
};

#endif /* PROFILEZOOM_H_ */
