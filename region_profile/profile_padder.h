/*
 * profile_padder.h
 *
 *  Created on: Jun 27, 2011
 *      Author: xin
 */

#ifndef PROFILE_PADDER_H_
#define PROFILE_PADDER_H_

#include <vector>
#include <stdint.h>

template<typename T>
class profile_padder {
public:
    static void pad(std::vector<T>& region) {

        uint32_t start_ind = 0;
        uint32_t end_ind = 0;
        T current_value;
        T pad_val;
        bool needTheEnd = false;
        int pad_region_length;

        for (uint32_t i = 0; i < region.size(); i++) {
            current_value = region[i];
            if (current_value != 0) {
                if (!needTheEnd) {
                    start_ind = i;
                    needTheEnd = true;
                } else {
                    end_ind = i;
                    if (end_ind > start_ind + 1) {

                        pad_val =
                                (T) ((region[start_ind] + region[end_ind]) / 2);

                        for (uint32_t ind = start_ind + 1; ind < end_ind;
                                ind++) {
                            region[ind] = pad_val;
                        }
                    }
                    needTheEnd = true;
                    start_ind = i;
                }
            }
        }
    }
};

#endif /* PROFILE_PADDER_H_ */
