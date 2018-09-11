/*
 * Strand.cpp
 *
 *  Created on: May 3, 2012
 *      Author: xin
 */

#include "Strand.h"

namespace reads {

Strand::Strand():mDir(true) {

}

Strand::~Strand() {

}

Strand::Strand(bool pos) :
        mDir(pos) {
}

bool Strand::isPos() const {
    return mDir ? true : false;
}

bool Strand::isNeg() const {
    return mDir ? false : true;
}

void Strand::toPos() {
    mDir = true;
}

void Strand::operator =(bool dir) {
    if (dir) {
        toPos();
    } else {
        toNeg();
    }
}

bool Strand::operator ==(const Strand& rhs) const {
	return mDir==rhs.mDir;
}

Strand::Strand(const Strand& str):mDir(str.mDir) {}

void Strand::toNeg() {
	mDir = false;
}

}

