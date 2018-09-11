/*
 * BamFileSortOrderAux.h
 *
 *  Created on: May 30, 2012
 *      Author: xfeng
 */

#ifndef BAMFILESORTORDERAUX_H_
#define BAMFILESORTORDERAUX_H_
#include "bamtools/BamAux.h"

namespace BamTools {
bool operator>(const BamTools::BamAlignment& bam1,
        const BamTools::BamAlignment& bam2);
bool operator<(const BamTools::BamAlignment& bam1,
        const BamTools::BamAlignment& bam2);
}

namespace parser {

bool readsSayCoordSorted(const std::string& file, uint32_t linesToTest);
bool readsSayReadNameSorted(const std::string& file, uint32_t linesToTest);
bool headerSaysCoordSorted(const std::string& file);
} /* namespace parser */
#endif /* BAMFILESORTORDERAUX_H_ */
