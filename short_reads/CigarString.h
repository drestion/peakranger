/*
 * CigarString.h
 *
 *  Created on: May 11, 2012
 *      Author: xin
 */

#ifndef CIGARSTRING_H_
#define CIGARSTRING_H_

#include "concepts/RegionInt32.h"
#include <stdint.h>
namespace reads {

class CigarString {
    friend std::ostream& operator<<(std::ostream& os, const CigarString& g);
public:
    CigarString();
    CigarString(const size_t len, const int32_t offset);
    virtual ~CigarString();

    virtual ranger::concepts::RegionInt32 toRegion() const;
    virtual std::string toString() const;

    size_t getLength() const;
    int32_t getOffset() const;

    void setLength(size_t mLength);
    void setOffset(int32_t mOffset);

protected:
    size_t mLength;
    int32_t mOffset;
};


class CigarS: public reads::CigarString {
public:
    CigarS();
    virtual ~CigarS();
    CigarS(const size_t& len, const size_t& offset);
    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};


class CigarP: public reads::CigarString {
public:
    CigarP();
    virtual ~CigarP();
    CigarP(const size_t& len, const size_t& offset);
    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};


class CigarN: public reads::CigarString {
public:
    CigarN();
    virtual ~CigarN();
    CigarN(const size_t& len, const size_t& offset);
    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};


class CigarM : public CigarString{
public:
    CigarM();
    CigarM(const size_t& len, const size_t& offset);
    virtual ~CigarM();

    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};


class CigarI: public reads::CigarString {
public:
    CigarI();
    virtual ~CigarI();
    CigarI(const size_t& len, const size_t& offset);
    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};


class CigarH: public reads::CigarString {
public:
    CigarH();
    virtual ~CigarH();
    CigarH(const size_t& len, const size_t& offset);
    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};


class CigarD: public reads::CigarString {
public:
    CigarD();
    virtual ~CigarD();
    CigarD(const size_t& len, const size_t& offset);
    ranger::concepts::RegionInt32 toRegion() const;
    std::string toString() const;
};
} /* namespace reads */
#endif /* CIGARSTRING_H_ */
