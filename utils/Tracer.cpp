/*
 * Tracer.cpp
 *
 *  Created on: Mar 23, 2012
 *      Author: xin
 */

#include "Tracer.h"
using namespace std;

namespace utils {

Tracer::Tracer() :
        os(cout), verbose(true) {

}

Tracer::Tracer(std::ostream & os, bool verb) :
        os(os), verbose(verb) {
}

Tracer::~Tracer() {

}

bool Tracer::isVerbose() const {
    return verbose;
}

void Tracer::setVerbose(bool verbose) {
    this->verbose = verbose;
}

}

utils::TimeStampTracer::TimeStampTracer(): os(cout), verbose(true) {
}

utils::TimeStampTracer::TimeStampTracer(std::ostream& os, bool verb) :
        os(os), verbose(verb) {
}



/* namespace utils */
