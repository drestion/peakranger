/*
 * BarCounter.h
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#ifndef BARCOUNTER_H_
#define BARCOUNTER_H_
#include <boost/icl/interval_map.hpp>

namespace boost {
namespace icl {

typedef interval_map<int32_t, int32_t> BarCounter;
typedef discrete_interval<int32_t> ICLInterval;
void operator+=(BarCounter& b1, const BarCounter& b2);

}
}
#endif /* BARCOUNTER_H_ */
