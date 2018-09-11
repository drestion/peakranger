/*
 * Region.h
 *
 *  Created on: May 10, 2012
 *      Author: xinfeng
 */

#ifndef REGIONCONCEPTS_H_
#define REGIONCONCEPTS_H_

#include "common/stl_header.h"
namespace ranger {
namespace concepts {

template<typename T>
class Region;

template<typename T>
std::ostream& operator<<(std::ostream& os, Region<T>& r) {
    os << r.getL() << "\t" << r.getR();
    return os;
}

template<typename LessThanComparable>
class Region {

public:

    Region();
    Region(const LessThanComparable& l, const LessThanComparable& r);

    std::ostream& operator<<(std::ostream& os);

    bool operator<(const Region<LessThanComparable>& rhs) const;
    bool operator>(const Region<LessThanComparable>& rhs) const;
    bool operator==(const Region<LessThanComparable>& rhs) const;

    const Region<LessThanComparable>& operator+(const int32_t& offset);
    bool operator[](const Region<LessThanComparable>& rhs) const;
    bool overlaps(const Region<LessThanComparable>& rhs) const;

    Region<LessThanComparable>& operator=(
            const Region<LessThanComparable>& rhs);

    LessThanComparable getL() const;
    LessThanComparable getR() const;
    void setL(const LessThanComparable& l);
    void setR(const LessThanComparable& r);

private:
    LessThanComparable mL;
    LessThanComparable mR;
};

template<typename LessThanComparable>
concepts::Region<LessThanComparable>::Region() :
        mL(), mR() {
}

template<typename LessThanComparable>
concepts::Region<LessThanComparable>::Region(const LessThanComparable& l,
        const LessThanComparable& r) :
        mL(l), mR(r) {
}

template<typename LessThanComparable>
bool concepts::Region<LessThanComparable>::operator <(
        const Region<LessThanComparable> & rhs) const {
    return mR < rhs.mL;
}

template<typename LessThanComparable>
bool concepts::Region<LessThanComparable>::operator >(
        const Region<LessThanComparable> & rhs) const {
    return mL > rhs.mR;
}

template<typename LessThanComparable>
bool concepts::Region<LessThanComparable>::operator ==(
        const Region<LessThanComparable> & rhs) const {
    return mL == rhs.mL && mR == rhs.mR;
}

template<typename LessThanComparable>
bool concepts::Region<LessThanComparable>::operator [](
        const Region<LessThanComparable> & rhs) const {
    return !(*this < rhs) && !(*this > rhs);
}

template<typename LessThanComparable>
bool concepts::Region<LessThanComparable>::overlaps(
        const Region<LessThanComparable> & rhs) const {
    return !(*this < rhs) && !(*this > rhs);
}

template<typename LessThanComparable>
Region<LessThanComparable> & concepts::Region<LessThanComparable>::operator =(
        const Region<LessThanComparable> & rhs) {
    mL = rhs.mL;
    mR = rhs.mR;
    return *this;
}

template<typename LessThanComparable>
LessThanComparable concepts::Region<LessThanComparable>::getL() const {
    return mL;
}

template<typename LessThanComparable>
void concepts::Region<LessThanComparable>::setL(const LessThanComparable& l) {
    mL = l;
}

template<typename LessThanComparable>
const Region<LessThanComparable> & Region<LessThanComparable>::operator +(
        const int32_t & offset) {
    mL += offset;
    mR += offset;
    return *this;
}

template<typename LessThanComparable>
std::ostream& Region<LessThanComparable>::operator <<(std::ostream& os) {
    os << getL() << "\t" << getR();
    return os;
}

template<typename LessThanComparable>
void concepts::Region<LessThanComparable>::setR(const LessThanComparable& r) {
    mR = r;
}

template<typename LessThanComparable>
LessThanComparable concepts::Region<LessThanComparable>::getR() const {
    return mR;
}

}/* namespace concepts */
}/* namespace ranger */
#endif /* REGIONCONCEPTS_H_ */

