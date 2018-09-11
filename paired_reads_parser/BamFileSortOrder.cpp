/*
 * BamFileSortOrder.cpp
 *
 *  Created on: May 30, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/BamFileSortOrder.h"
#include "bamtools/BamAux.h"
#include "bamtools/BamReader.h"
#include "utils/Guarded.h"
#include "common/boost_header.h"
#include "BamFileSortOrderAux.h"
using namespace BamTools;
using namespace utils;
using namespace std;
namespace parser {

BamFileSortOrder::BamFileSortOrder() {

}

BamFileSortOrder::~BamFileSortOrder() {

}

bool BamFileSortOrder::match(const std::string& file, uint32_t n) {

    return true;
}

bool BamFileCoordinateSortedOrder::match(const std::string& file, uint32_t n) {

//    return headerSaysCoordSorted(file) || readsSayCoordSorted(file, n);
    return  readsSayCoordSorted(file, n);
}

bool BamFileCoordinateUnSortedOrder::match(const std::string& file,
        uint32_t n) {
    return !readsSayCoordSorted(file, n);
}

bool BamFileReadNameSortedOrder::match(const std::string& file, uint32_t n) {
    return readsSayReadNameSorted(file, n);
}

bool BamFileReadNameUnSortedOrder::match(const std::string& file, uint32_t n) {
    return !readsSayReadNameSorted(file, n);
}

}

/* namespace parser */
