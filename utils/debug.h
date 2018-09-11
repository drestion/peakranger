/*
 * debug.h
 *
 *  Created on: Mar 22, 2012
 *      Author: xin
 */

#ifndef DEBUG_H_
#define DEBUG_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/foreach.hpp>
#include <stdint.h>
#define foreach BOOST_FOREACH

#define COUT(msg,obj) \
    std::cout <<msg<<"\t"<<obj<<"\n";
#define QUIT(msg) \
    std::cout << msg <<std::endl;exit(0)

namespace ranger_debug {
template<typename T>
void dumpArray(T* arr,
               size_t size,
               const char* filename,
               const char* delim = "\n") {
    std::ofstream os(filename);
    for (size_t i = 0; i < size; i++) {
        os << i << "\t" << arr[i] << delim;
    }
}
template<typename T>
void dumpArray(const T& arr,
               const char* filename,
               const char* delim = "\n") {
    std::ofstream os(filename);
    typename T::const_iterator it = arr.begin();
    for (size_t i = 0; it != arr.end(); ++it) {
        os << i++ << "\t" << *it << delim;
    }
}

template<typename T, typename M>
void dumpArray(const T& arr,
               const char* filename,
               const char* delim = "\n") {
    std::ofstream os(filename);
    typename T::const_iterator it = arr.begin();
    for (size_t i = 0; it != arr.end(); ++it) {
        os << i++ << "\t" << M::map(*it) << delim;
    }
}

template<typename T, typename M>
void dumpArray(T* arr,
               size_t size,
               const char* filename,
               const char* delim = "\n") {
    std::ofstream os(filename);
    for (size_t i = 0; i < size; i++) {
        os << i << "\t" << M::map(arr[i]) << delim;
    }
}
}

#endif /* DEBUG_H_ */
