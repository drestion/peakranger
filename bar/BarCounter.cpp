/*
 * BarCounter.cpp
 *
 *  Created on: Jun 20, 2012
 *      Author: xfeng
 */

#include "BarCounter.h"
#include <iostream>
namespace boost {
namespace icl {

void operator+=(BarCounter& b1, const BarCounter& b2) {
    BarCounter::const_iterator it;
    if (&b1 == &b2) {
//        //will this comparison work?
        BarCounter bb2(b2);
        for (it = bb2.begin(); it != bb2.end(); it++) {
            b1 += make_pair((*it).first, (*it).second);
        }
    } else {
        for (it = b2.begin(); it != b2.end(); it++) {
            b1 += make_pair((*it).first, (*it).second);
        }
    }
}

}
}
