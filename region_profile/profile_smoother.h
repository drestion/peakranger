/*
 * profile_smoother.h
 *
 *  Created on: Jun 23, 2011
 *      Author: xin
 */

#ifndef PROFILE_SMOOTHER_H_
#define PROFILE_SMOOTHER_H_

#include <string>
#include <iostream>
#include <algorithm>
#include <stdint.h>

#include "utils/assert_helpers.h"
#include "utils/exceptions.h"

template<typename T>
class profile_smoother {

private:

    static inline void _average(std::vector<T>& vec, std::vector<double>& vout,
            size_t bandwidth) {
        uint32_t halfspan = bandwidth / 2;
        if (vec.size() < halfspan * 2)
            throw RangerException("In _average, vec.size < halfspan*2");

        vout.resize(vec.size() - halfspan * 2);

        for (size_t i = 0; i < vout.size(); i++) {
            for (int32_t j = -halfspan; j < (int32_t) (halfspan + 1); j++) {

                vout.at(i) += vec.at(j + i + halfspan);
            }

            vout[i] /= 1.0 * (halfspan * 2 + 1);
        }
    }

    static inline void _arrayPlus(std::vector<T>& array, uint32_t start,
            uint32_t end) {
        assert(end<=array.size());
        assert_lt(start, end)
        for (uint32_t i = start; i < end; i++) {
            array[i] += 1;
        }
    }

    static inline void _arrayDevideForOddMember(std::vector<double>& array,
            std::vector<double>& result, uint32_t size) {
        uint32_t mid = size;
        if ((size % 2))
            mid = (size + 1) / 2;
        else
            mid = size / 2;
        result.resize(mid);
        int32_t k = -2;
        for (uint32_t i = 0; i < mid; i++) {
            k += 2;
            result[i] = array.at(k) * 1.0 / (k + 1);
        }
    }

    static inline void _cumsum(std::vector<T>& vec, std::vector<double>& vout,
            uint32_t size) {
        vout.resize(size);
        assert_lt(size-1, vec.size())
        double cums = 0;
        for (uint32_t i = 0; i < size; i++) {
            vout[i] = cums + vec[i];
            cums += vec[i];
        }
    }

    static inline void _reverse_cumsum(std::vector<T>& data,
            std::vector<double>& vout, uint32_t offset, uint32_t size) {

        vout.resize(size);
        double cums = 0;
        uint32_t ind = size - 1;
        assert_lt(ind+offset, data.size())
        for (uint32_t i = 0; i < size; i++) {
            vout[i] = cums + data[offset + ind];
            cums += data[offset + ind];
            ind--;
        }
    }

public:
    static void smooth(std::vector<T>& profile, uint32_t bandwidth) {
        if (bandwidth <= 1) {
            //No need to smooth
            return;
        }

        uint32_t length = profile.size();
        bandwidth = bandwidth > length ? length : bandwidth;
        bandwidth = bandwidth - 1 + (bandwidth % 2); // force it to be odd;
        std::vector<double> middle, cbegin, cbegin_h, cend, cend_h;

        _average(profile, middle, bandwidth);

        _cumsum(profile, cbegin, bandwidth - 2);

        _arrayDevideForOddMember(cbegin, cbegin_h, cbegin.size());

        _reverse_cumsum(profile, cend, length - bandwidth + 2, bandwidth - 2);

        _arrayDevideForOddMember(cend, cend_h, cend.size());

        std::reverse(cend_h.begin(), cend_h.end());
        assert_lt(cbegin_h.size()-1, profile.size())
        for (size_t i = 0; i < cbegin_h.size(); i++) {
            profile[i] = cbegin_h[i];
        }
        assert_lt(cbegin_h.size()+middle.size()-1, profile.size())
        for (size_t i = 0; i < middle.size(); i++) {
            profile[i + cbegin_h.size()] = middle[i];
        }
        assert_lt(cbegin_h.size()+middle.size()+cend_h.size()-1, profile.size())
        for (size_t i = 0; i < cend_h.size(); i++) {
            profile[i + cbegin_h.size() + middle.size()] = cend_h[i];
        }

    }
};
#endif /* PROFILE_SMOOTHER_H_ */
