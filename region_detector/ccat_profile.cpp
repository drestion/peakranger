/*
 * ccat_profile.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#include "ccat_profile.h"
#include "ccat_profile_helper.h"
#include "common/boost_header.h"

#include <functional>
#include <boost/lambda/lambda.hpp>
#include <boost/bind.hpp>
#include "math/histogram.h"
#include "utils/stl_helper.h"
using namespace std;
using namespace boost::lambda;
using namespace boost;

namespace ccat_aux {

ccat_profile_t::ccat_profile_t() :
        fragmentSize(1), movingStep(1), profileLength(1), chromSize(1) {

}

ccat_profile_t::~ccat_profile_t() {

}

void ccat_profile_t::ResampleProfile(const std::vector<size_t> & vals, std::vector<size_t> & result,
        size_t profileLength, double ratio, bool isNegStrand) {
    EnforceMaxIndex(profileLength);
    assert(result.size() == profileLength);

    size_t ind;

    foreach(size_t val, vals) {

        if (isNegStrand) {
            ind = NegReadMapper(val);
        } else {
            ind = PosReadMapper(val);
        }
        if (!ccat_aux::rsRatio(ratio)) {

            result[ind] += 1;

        }

    }
}

void ccat_profile_t::rsPosReadMapper(size_t read, double ratio, size_t& mapped) {

}

void ccat_profile_t::rsNegReadMapper(size_t read, double ratio, size_t& mapped) {
    //size_t res = NegReadMapper(read);
    NegReadMapper(read);

}

void ccat_profile_t::EnforceMaxIndex(size_t profileLength) {
    this->profileLength = profileLength;
}

void ccat_profile_t::GetProfile(const std::vector<size_t> & vals, std::vector<size_t> & result, size_t profileLength,
        bool isNegStrand) {

    EnforceMaxIndex(profileLength);
    assert_gt(profileLength, 2);
    if (result.size() != profileLength) {
        ranger_stl_aux::resetVecToSize<size_t>(result, 0, profileLength);
    }
    if (isNegStrand) {
        ranger_math::getHist(vals, std::bind1st(std::mem_fun(&ccat_profile_t::NegReadMapper), this), result);
    } else {
        ranger_math::getHist(vals, std::bind1st(std::mem_fun(&ccat_profile_t::PosReadMapper), this), result);
    }
}

size_t ccat_profile_t::PosReadMapper(size_t read) {

    size_t tmpIndex = (read + fragmentSize / 2) / movingStep;

    tmpIndex = tmpIndex >= profileLength ? profileLength - 1 : tmpIndex;
    return tmpIndex;
}

size_t ccat_profile_t::NegReadMapper(size_t read) {

    size_t tmpIndex;
    if (read < fragmentSize / 2) {
        tmpIndex = 0;
    } else {
        tmpIndex = (read - fragmentSize / 2) / movingStep;
    }

    tmpIndex = tmpIndex >= profileLength ? profileLength - 1 : tmpIndex;
    return tmpIndex;
}

size_t ccat_profile_t::getFragmentSize() const {
    return fragmentSize;
}

size_t ccat_profile_t::getMovingStep() const {
    return movingStep;
}

size_t ccat_profile_t::getProfileLength() const {
    return profileLength;
}

void ccat_profile_t::setFragmentSize(size_t fragmentSize) {
    this->fragmentSize = fragmentSize;
}

void ccat_profile_t::setMovingStep(size_t ms) {
    this->movingStep = ms;
}

void ccat_profile_t::setProfileLength(size_t profileLength) {
    this->profileLength = profileLength;
}

size_t ccat_profile_t::getChromSize() const {
    return chromSize;
}

void ccat_profile_t::setChromSize(size_t chromSize) {
    this->chromSize = chromSize;
}

/* namespace ccat_aux */
}

