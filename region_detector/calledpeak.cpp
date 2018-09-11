/*
 * calledpeak.cpp
 *
 *  Created on: Dec 9, 2011
 *      Author: xfeng
 */

#include "calledpeak.h"

std::ostream & operator <<(std::ostream & os, const called_peak & pk) {
    os << pk.first << "\t" << pk.second << "\t" << pk.p << "\t" << pk.q << "\t"
            <<

            pk.treads << "\t" << pk.creads << "\t";
    return os;
}

