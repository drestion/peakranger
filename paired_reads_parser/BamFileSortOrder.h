/*
 * BamFileSortOrder.h
 *
 *  Created on: May 30, 2012
 *      Author: xfeng
 */

#ifndef BAMFILESORTORDER_H_
#define BAMFILESORTORDER_H_

#include <string>
#include <stdint.h>
namespace parser {

enum BamFileSortOrderType {
    coordinate_sorted = 0x1,
    readname_sorted = 0x2,
    coordinate_unsorted = 0x4,
    readname_unsorted = 0x8,
    unknown_order = 0x0
};

/*
 * Determine the sort order of bam
 * file by checking required lines
 */
class BamFileSortOrder {
public:
    BamFileSortOrder();
    virtual ~BamFileSortOrder();

    virtual bool match(const std::string& file, uint32_t linesToTest);
};

class BamFileCoordinateSortedOrder: public BamFileSortOrder {
public:
    bool match(const std::string& file, uint32_t linesToTest);
};

class BamFileCoordinateUnSortedOrder: public BamFileSortOrder {
public:
    bool match(const std::string& file, uint32_t linesToTest);
};

class BamFileReadNameSortedOrder: public BamFileSortOrder {
public:
    bool match(const std::string& file, uint32_t linesToTest);
};

class BamFileReadNameUnSortedOrder: public BamFileSortOrder {
public:
    bool match(const std::string& file, uint32_t linesToTest);
};

} /* namespace parser */
#endif /* BAMFILESORTORDER_H_ */
