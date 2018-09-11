/*
 * BamFileSortOrderDetector.h
 *
 *  Created on: Jun 5, 2012
 *      Author: xfeng
 */

#ifndef BAMFILESORTORDERDETECTOR_H_
#define BAMFILESORTORDERDETECTOR_H_
#include "paired_reads_parser/BamFileSortOrder.h"

namespace parser {

class BamFileSortOrderDetector {
public:
    BamFileSortOrderDetector();
    virtual ~BamFileSortOrderDetector();
    uint32_t getType(const std::string& file, uint32_t linesToUse);

};

} /* namespace parser */
#endif /* BAMFILESORTORDERDETECTOR_H_ */
