/*
 * wig.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: xfeng
 */

#include "wig.h"


uint32_t wig::getP() const
{
    return _p;
}

double wig::getS() const
{
    return _s;
}

void wig::setP(uint32_t p)
{
    _p = p;
}

void wig::setS(double s)
{
    _s = s;
}

wig::wig()
{
}

wig::wig(uint32_t p, double s)
{
    _p = p;
    _s = s;
}

bool wig::operator<(const wig & w)
{
    return this->_p < w._p;
}

bool wig::compa(const wig & w1, const wig & w2)
{
    return w1._p < w2._p;
}





