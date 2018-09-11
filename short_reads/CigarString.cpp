/*
 * CigarString.cpp
 *
 *  Created on: May 11, 2012
 *      Author: xin
 */

#include "CigarString.h"
using namespace ranger::concepts;

namespace reads {

CigarString::CigarString():mLength(0),mOffset(0) {

}

CigarString::CigarString(const size_t len, const int32_t offset) :
        mLength(len), mOffset(offset) {
}

CigarString::~CigarString() {

}

size_t CigarString::getLength() const {
    return mLength;
}

RegionInt32 CigarString::toRegion() const {
    return RegionInt32(getOffset(), getOffset() + getLength());
}

void CigarString::setLength(size_t mLength) {
    this->mLength = mLength;
}

int32_t CigarString::getOffset() const {
    return mOffset;
}

std::ostream& operator <<(std::ostream& os, const CigarString& g) {
    os << g.toString();
    return os;
}

void CigarString::setOffset(int32_t mOffset) {
    this->mOffset = mOffset;
}

std::string CigarString::toString() const {
    return "";
}

CigarS::CigarS() {


}

CigarS::~CigarS() {

}

ranger::concepts::RegionInt32 CigarS::toRegion() const {
    return RegionInt32(0, 0);
}

CigarS::CigarS(const size_t& len, const size_t& offset) {
    setLength(len);
    setOffset(offset);
}

std::string CigarS::toString() const {
    std::stringstream ss;
    ss << getLength() << 'S';
    return ss.str();
}

CigarP::CigarP() {

}

CigarP::~CigarP() {

}

ranger::concepts::RegionInt32 CigarP::toRegion() const {
    return RegionInt32(0, 0);
}

CigarP::CigarP(const size_t& len, const size_t& offset) {
    setLength(len);
    setOffset(offset);
}

std::string CigarP::toString() const {
    std::stringstream ss;
    ss << getLength() << 'P';
    return ss.str();
}


CigarN::CigarN() {


}

CigarN::~CigarN() {

}

ranger::concepts::RegionInt32 CigarN::toRegion() const {
    return RegionInt32(0, 0);
}

CigarN::CigarN(const size_t& len, const size_t& offset) {
    setLength(len);
    setOffset(offset);
}

std::string CigarN::toString() const {
    std::stringstream ss;
    ss << getLength() << 'N';
    return ss.str();
}


CigarM::CigarM() {

}

CigarM::~CigarM() {

}

CigarM::CigarM(const size_t & strLen, const size_t& offset) {
    mOffset = offset;
    mLength = strLen;
}

RegionInt32 CigarM::toRegion() const {
    //0-based
    return RegionInt32(mOffset, mOffset + mLength);
}

std::string CigarM::toString() const {
    std::stringstream ss;
    ss << getLength() << 'M';
    return ss.str();
}


CigarI::CigarI() {


}

CigarI::~CigarI() {

}

ranger::concepts::RegionInt32 CigarI::toRegion() const {
    return RegionInt32(0, 0);
}

CigarI::CigarI(const size_t& len, const size_t& offset) {
    setLength(len);
    setOffset(offset);
}

std::string CigarI::toString() const {
    std::stringstream ss;
    ss << getLength() << 'I';
    return ss.str();
}


CigarH::CigarH() {


}

CigarH::~CigarH() {

}

ranger::concepts::RegionInt32 CigarH::toRegion() const {
    return RegionInt32(0, 0);
}

CigarH::CigarH(const size_t& len, const size_t& offset) {
    setLength(len);
    setOffset(offset);
}

std::string CigarH::toString() const {
    std::stringstream ss;
    ss << getLength() << 'H';
    return ss.str();
}


CigarD::CigarD() {
}

CigarD::~CigarD() {
}

ranger::concepts::RegionInt32 CigarD::toRegion() const {
    return RegionInt32(0, 0);
}

CigarD::CigarD(const size_t& len, const size_t& offset) {
    setLength(len);
    setOffset(offset);
}

std::string CigarD::toString() const {
    std::stringstream ss;
    ss << getLength() << 'D';
    return ss.str();
}

} /* namespace reads */
