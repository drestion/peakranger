/*
 * ReadPair.h
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#ifndef READPAIR_H_
#define READPAIR_H_

#include "common/stl_header.h"

namespace reads {

template<typename R>
class ReadPair;

template<typename R>
std::ostream& operator<<(std::ostream& os, ReadPair<R>& rhs) {
    os << rhs.r1();
    os << rhs.r2();
    return os;
}

template<typename R>
class ReadPair {
public:
    ReadPair();
    ReadPair(const R& r1, const R& r2);
    virtual ~ReadPair();

    R& r1();
    R& r2();
    const R& r1() const;
    const R& r2() const;
    void add(const R& r1, const R& r2);

    bool operator==(const ReadPair<R>& rhs) const {
        return mR1 == rhs.mR1 && mR2 == mR2;
    }

private:
    R mR1;
    R mR2;
};

}

template<typename R> reads::ReadPair<R>::ReadPair() :
        mR1(), mR2() {
}

template<typename R> reads::ReadPair<R>::ReadPair(const R & r1, const R & r2) :
        mR1(r1), mR2(r2) {
}

template<typename R> reads::ReadPair<R>::~ReadPair() {
}

template<typename R> R & reads::ReadPair<R>::r1() {
    return mR1;
}

template<typename R> R & reads::ReadPair<R>::r2() {
    return mR2;
}

template<typename R> void reads::ReadPair<R>::add(const R & r1, const R & r2) {
    mR1 = r1;
    mR2 = r2;
}

template<typename R>
inline const R& reads::ReadPair<R>::r1() const {
    return mR1;
}

template<typename R>
inline const R& reads::ReadPair<R>::r2() const {
    return mR2;
}


/* namespace reads */
#endif /* READPAIR_H_ */
