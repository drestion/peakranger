/*
 * histogram.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <vector>
#include <stdint.h>
#include <algorithm>

namespace ranger_math {

template<typename C, typename T, typename M>
void getHist(const T& vals,
             M pf,
             C& res) {

    typename T::const_iterator it;
    for (it = vals.begin(); it != vals.end(); it++) {
        res.at(pf(*it))++;
        }

    }

template<typename C,typename T,typename M,typename R>
void getHist(const T& vals, C& res) {

    typename T::const_iterator it;
    for (it = vals.begin(); it != vals.end(); it++) {
        R::map(res.at(M::map(*it)))++;
    }

}

} /* namespace ranger_math */
#endif /* HISTOGRAM_H_ */
