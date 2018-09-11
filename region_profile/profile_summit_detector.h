/*
 * detector.h
 *
 *  Created on: Jun 22, 2011
 *      Author: xin
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <limits.h>
#include <math.h>

#include "utils/assert_helpers.h"
#include "utils/exceptions.h"

#define MAX_NUM_OF_SUBPEAKS 100000

class profile_summit_detector {
public:
    template<typename T>
    static void detect_sub_peaks(std::vector<T> smoothed_profile, double delta,
            std::vector<uint32_t>& max_pos, std::vector<uint32_t>& min_pos) {
        {
            rt_assert_gt(delta, 0)
            rt_assert_lt(delta, 1.0000000001)
            double mx = -INFINITY;
            double mn = INFINITY;
            double threshold = (double) (*max_element(smoothed_profile.begin(),
                    smoothed_profile.end()));
            threshold *= delta;
            int look_for_max = 1;
            uint32_t mxpos = 0;
            uint32_t mnpos = 0;
            uint32_t kk = 0;
            uint32_t jj = 0;
            uint32_t profile_length = smoothed_profile.size();

            for (uint32_t i = 0; i < profile_length; i++) {
                if (kk >= MAX_NUM_OF_SUBPEAKS || jj >= MAX_NUM_OF_SUBPEAKS) {
                    throw RangerException(
                            "in summit detector: Too many subpeaks were detected.");
                }
                if (smoothed_profile[i] > mx) {
                    mx = smoothed_profile[i];
                    mxpos = i;
                }
                if (smoothed_profile[i] < mn) {
                    mn = smoothed_profile[i];
                    mnpos = i;
                }
                if (look_for_max) {
                    if (smoothed_profile[i] < (mx - threshold)) {
                        max_pos.push_back(mxpos);
                        jj++;
                        mn = smoothed_profile[i];
                        mnpos = i;
                        look_for_max = 0;
                    }
                } else {
                    if (smoothed_profile[i] > (mn + threshold)) {
                        min_pos.push_back(mnpos);
                        kk++;
                        mx = smoothed_profile[i];
                        mxpos = i;
                        look_for_max = 1;
                    }
                }
            }

            if (jj < 1) {
                max_pos.push_back(mxpos);
                min_pos.push_back(mnpos);
            }
            if (kk < 1) {
                min_pos.push_back(mnpos);
            }
        }
    }
};
#endif /* DETECTOR_H_ */
