/*
 * ccat_profile.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef CCAT_PROFILE_H_
#define CCAT_PROFILE_H_

#include "common/stl_header.h"
#include "common/ranger_debug.h"
namespace ccat_aux {

class ccat_profile_t {
public:
    ccat_profile_t();
    virtual ~ccat_profile_t();

    template<typename T>
    void SmoothProfile(std::vector<T> & array, size_t len, size_t halfWinSize) {
        //SmoothProfile: Smooth a profile by summing over a sliding window

        size_t i;
        int sum;
        std::vector<T> tmp(len, 0);
        assert(len > halfWinSize);
        sum = 0;
        for (i = 0; i < halfWinSize * 2 + 1; i++) {
            sum += array[i];
        }

        for (i = 0; i <= halfWinSize; i++) {
            tmp[i] = sum;
        }

        for (i = halfWinSize + 1; i < len - halfWinSize; i++) {
            tmp[i] = tmp[i - 1] - array[i - halfWinSize - 1] + array[i + halfWinSize];
        }

        for (i = len - halfWinSize; i < len; i++) {
            tmp[i] = tmp[i - 1];
        }

        array.swap(tmp);
    }

    void GetProfile(const std::vector<size_t>& tags, std::vector<size_t>& result, size_t profileLength,
            bool isNegStrand = false);
    void ResampleProfile(const std::vector<size_t>& tags, std::vector<size_t>& result, size_t profileLength,
            double ratio, bool isNegStrand = false);
    size_t getFragmentSize() const;
    size_t getMovingStep() const;
    size_t getProfileLength() const;
    void setFragmentSize(size_t fragmentSize);
    void setMovingStep(size_t movingStep);
    void setProfileLength(size_t profileLength);
    size_t getChromSize() const;
    void setChromSize(size_t chromSize);

private:
    size_t PosReadMapper(size_t read);
    size_t NegReadMapper(size_t read);
    void rsPosReadMapper(size_t read, double ratio, size_t& mapped);
    void rsNegReadMapper(size_t read, double ratio, size_t& mapped);
    void EnforceMaxIndex(size_t profileLength);
    size_t fragmentSize;
    size_t movingStep;
    size_t profileLength;
    size_t chromSize;
};

} /* namespace ccat_aux */
#endif /* CCAT_PROFILE_H_ */
