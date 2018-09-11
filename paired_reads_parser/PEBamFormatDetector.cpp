/*
 * PEBamFormatDetector.cpp
 *
 *  Created on: Jun 6, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/PEBamFormatDetector.h"
#include "paired_reads_parser/BamParserAux.h"
#include "short_reads/bamFormatAux.h"
#include "short_reads/BlockedBamFormatAux.h"
using namespace BamTools;
using namespace std;
namespace parser {

PEBamFormatDetector::PEBamFormatDetector() {

}

PEBamFormatDetector::~PEBamFormatDetector() {

}

bool PEBamFormatDetector::isPairEndBamFile(const std::string& file,
        uint32_t linesToTest) {
    vector<BamAlignment> bamsToCheck;
    fetchLines(bamsToCheck, linesToTest, file);
//    vector<BamAlignment>::iterator it;
//    it = std::find_if(bamsToCheck.begin(), bamsToCheck.end(), reads::isSERead);
//    if (it != bamsToCheck.end()) {
//        cout << "[PEBamFormatDetector::isPairEndBamFile]" << *it << "\n";
//    }
    return std::find_if(bamsToCheck.begin(), bamsToCheck.end(), reads::isSERead)
            == bamsToCheck.end();
}

} /* namespace parser */
