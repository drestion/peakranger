/*
 * stl_helper.h
 *
 *  Created on: Mar 14, 2012
 *      Author: xfeng
 */

#ifndef STL_HELPER_H_
#define STL_HELPER_H_

#include <iostream>
#include <functional>
#include <vector>
namespace ranger_stl_aux {

template<class T> struct print_sep_empty:public std::unary_function<T, void>
{
    print_sep_empty(std::ostream& out)
    : os(out), count(0) {
    }
    void operator()(T x) {
        os << x << ' ';
        ++count;
    }
    std::ostream& os;
    int count;
};

template<class T> struct print_sep_tab:public std::unary_function<T, void>
{
    print_sep_tab(std::ostream& out)
    : os(out), count(0) {
    }
    void operator()(T x) {
        os << x << '\t';
        ++count;
    }
    std::ostream& os;
    int count;
};

template<class T> struct print_sep_comma:public std::unary_function<T, void>
{
    print_sep_comma(std::ostream& out)
    : os(out), count(0) {
    }
    void operator()(T x) {
        os << x << '\t';
        ++count;
    }
    std::ostream& os;
    int count;
};

template<class T> void resetVec(std::vector<T>& vec, const T& defVal){
    std::vector<T> tmp(vec.size(),defVal);
    vec.swap(tmp);
}


template<class T> void resetVecToSize(std::vector<T>& vec, const T& defVal, size_t newSize=0){
    std::vector<T> tmp(newSize,defVal);
    vec.swap(tmp);
}

}

#endif /* STL_HELPER_H_ */
