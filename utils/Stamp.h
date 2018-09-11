/*
 * Stamp.h
 *
 *  Created on: Jul 5, 2012
 *      Author: tania
 */

#ifndef STAMP_H_
#define STAMP_H_
#include "utils/timer.h"
#include "utils/util_print.h"
#include <iostream>
namespace utils {

class Stamp {
public:

	static void citationAndDate(std::ostream& os);
    static void citationRangerCCATAndDate(std::ostream& os);
    static void citationRangerBCPAndDate(std::ostream& os);
};

} /* namespace utils */
#endif /* STAMP_H_ */
