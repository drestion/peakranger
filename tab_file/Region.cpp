/*
 * Region.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: xinfeng
 */

#include "Region.h"
#include <assert.h>


namespace tab_file {

Region::Region()
: l(0), r(0) {

}

Region::Region(loc_t l,
               loc_t r)
: l(l), r(r) {
    assert(l <= r);
}

Region::Region(const Region& r)
: l(r.l), r(r.r) {
}

Region::~Region() {

}

bool Region::operator <(const Region& rhs) const {
    return r < rhs.l;
}

bool Region::operator >(const Region& rhs) const {
    return l > rhs.r;
}

size_t Region::getL() const {
    return l;
}

void Region::setL(size_t l) {
    this->l = l;
}

std::ostream & operator <<(std::ostream & os,
                           const Region & tg)
                           {
    os << tg.l << ", " << tg.r;
    return os;
}

size_t Region::getR() const {
    return r;
}

bool Region::operator [](const Region& rhs) const {
    return rhs.l >= l && rhs.r <= r;
}

void Region::setR(loc_t r) {
    this->r = r;
}

bool Region::operator ==(const Region& rhs) const {
    return l == rhs.l && r == rhs.r;
}

Region& Region::operator =(const Region& rhs)
                           {
    Region tmp(rhs);
    swap(tmp);
    return *this;
}

bool Region::overlaps(const Region & rhs) const
                      {
    return !(*this < rhs) && !(*this > rhs);
}

void Region::swap(Region & rr)
                  {
    loc_t lll = l, rrr = r;
    l = rr.l;
    r = rr.r;
    rr.l = lll;
    rr.r = rrr;
}

/* namespace tab_file */
}

