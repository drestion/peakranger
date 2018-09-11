/*
 * Tracer.h
 *
 *  Created on: Mar 23, 2012
 *      Author: xin
 */

#ifndef TRACER_H_
#define TRACER_H_
#include "utils/timer.h"
#include <iostream>

namespace utils {
class TimeStampTracer {
public:
    TimeStampTracer();
    virtual ~TimeStampTracer(){}
    TimeStampTracer(std::ostream& os, bool verb = true);


    template<typename T>
    TimeStampTracer& operator<<(const T& msg);

    bool isVerbose() const;
    void setVerbose(bool verbose);
    void flush();
    template<typename T>
    static void trace_if(const T& msg, bool cond, std::ostream& os = std::cout);

private:
    std::ostream& os;
    bool verbose;
};

class Tracer {
public:
    Tracer();
    Tracer(std::ostream& os, bool verb = true);
    virtual ~Tracer();

    template<typename T>
    Tracer& operator<<(const T& msg);

    bool isVerbose() const;
    void setVerbose(bool verbose);

    template<typename T>
    static void trace_if(const T& msg, bool cond, std::ostream& os = std::cout);

private:
    std::ostream& os;
    bool verbose;
};

template<typename T>
Tracer& Tracer::operator<<(const T& msg) {
    if (verbose) {
        os << msg;

    }
    return *this;
}
template<typename T>
void Tracer::trace_if(const T& msg, bool cond, std::ostream& os) {
    if (cond)
        os << msg;
}

template<typename T>
TimeStampTracer& TimeStampTracer::operator<<(const T& msg) {
    if (verbose) {
//        logDate(os, false);
        //todo: should use std::endl to signal data ouput
        os  << msg;
        os.flush();
    }
    return *this;
}

}

/* namespace utils */
#endif /* TRACER_H_ */
