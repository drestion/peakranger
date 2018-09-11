/*
 * reads_aux.h
 *
 *  Created on: Apr 12, 2012
 *      Author: xin
 */

#ifndef READS_AUX_H_
#define READS_AUX_H_
#include <algorithm>
#include <utility>
namespace reads_aux {

template<typename Iter, typename B>
void reads_in(Iter lower, Iter higher, const B& binstart, const B& binend, Iter& res_l, Iter& res_h) {
    res_l = std::lower_bound(lower, higher, binstart);
    res_h = std::upper_bound(res_l, higher, binend);
}
//
//template<typename R>
//bool inRegion(const R)
//
//template<typename Iter,>
//void reads_in2(Iter lower, Iter higher, const B& region, Iter& res_l, Iter& res_h) {
//
//    std::pair<Iter,Iter> res = std::equal_range()
//    res_l = std::lower_bound(lower, higher, binstart);
//    res_h = std::upper_bound(res_l, higher, binend);
//}

}

#endif /* READS_AUX_H_ */
