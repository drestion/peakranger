/*
 * BamFileSortOrderDetector.cpp
 *
 *  Created on: Jun 5, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/BamFileSortOrderDetector.h"
#include "common/boost_header.h"
using namespace boost;
namespace parser {

BamFileSortOrderDetector::BamFileSortOrderDetector() {

}

BamFileSortOrderDetector::~BamFileSortOrderDetector() {

}

uint32_t BamFileSortOrderDetector::getType(const std::string& file,
        uint32_t ln) {
    uint32_t res(unknown_order);
    if (make_shared<BamFileCoordinateSortedOrder>()->match(file, ln)) {
        res |= coordinate_sorted;
    }
    if (make_shared<BamFileReadNameSortedOrder>()->match(file, ln)) {
        res |= readname_sorted;
    }
    if (make_shared<BamFileCoordinateUnSortedOrder>()->match(file, ln)) {
        res |= coordinate_unsorted;
    }
    if (make_shared<BamFileReadNameUnSortedOrder>()->match(file, ln)) {
        res |= readname_unsorted;
    }
    return res;
}

} /* namespace parser */
