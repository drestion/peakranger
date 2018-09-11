/*
 * Stamp.cpp
 *
 *  Created on: Jul 5, 2012
 *      Author: tania
 */

#include "utils/Stamp.h"

namespace utils {

void Stamp::citationAndDate(std::ostream& os) {
    utilprint::citation ct;
    ct.print_msg(os);
    logDate(os);
}

void Stamp::citationRangerCCATAndDate(std::ostream& os) {
    utilprint::citation ct;
    utilprint::ccatCitation ct2;
    ct.print_msg(os);
    ct2.print_msg(os);
    logDate(os);
}

} /* namespace utils */

void utils::Stamp::citationRangerBCPAndDate(std::ostream& os) {
	utilprint::citation ct;
	    utilprint::bcpCitation ct2;
	    ct.print_msg(os);
	    ct2.print_msg(os);
	    logDate(os);
}
