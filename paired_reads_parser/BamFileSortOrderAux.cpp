/*
 * BamFileSortOrderAux.cpp
 *
 *  Created on: May 30, 2012
 *      Author: xfeng
 */

#include "paired_reads_parser/BamFileSortOrderAux.h"
#include "paired_reads_parser/BamFileSortOrder.h"
#include "bamtools/BamReader.h"
#include "utils/Guarded.h"
#include "common/boost_header.h"
#include "short_reads/bamFormatAux.h"
#include "short_reads/BlockedBamFormatAux.h"
#include "paired_reads_parser/BamParserAux.h"
using namespace BamTools;
using namespace utils;
using namespace std;
namespace BamTools {
bool operator>(const BamTools::BamAlignment& bam1,
        const BamTools::BamAlignment& bam2) {
    return bam1.Position > bam2.Position;
}
bool operator <(const BamTools::BamAlignment& bam1,
        const BamTools::BamAlignment& bam2) {
    return bam1.Position < bam2.Position;
}

}
namespace parser {

bool readsSayCoordSorted(const string& file, uint32_t n) {
    assert(n > 1);
    vector<BamAlignment> bamsToCheck;
    fetchLines(bamsToCheck, n, file);

//    vector<BamAlignment>::iterator it;
//    it = std::adjacent_find(bamsToCheck.begin(), bamsToCheck.end(),
//            std::greater<BamAlignment>());
//    if (it != bamsToCheck.end()) {
//        cout << "[readsSayCoordSorted]" << *it << "\n";
//    }
    return std::adjacent_find(bamsToCheck.begin(), bamsToCheck.end(),
            std::greater<BamAlignment>()) == bamsToCheck.end();
}

bool readsSayReadNameSorted(const std::string& file, uint32_t n) {
    assert(n > 1);
    n = (n % 2 == 0) ? n : (n - 1);
    vector<BamAlignment> bamsToCheck;
    fetchLines(bamsToCheck, n, file);
    return std::adjacent_find(bamsToCheck.begin(), bamsToCheck.end(),
            reads::isDifferentReadName) == bamsToCheck.end();
}

bool headerSaysCoordSorted(const std::string& file) {
    BamReader bam;
    BamAlignment read;
    Guarded<FileNotGood> g(!(bam.Open(file)), file.c_str());
    const RefVector refvec = bam.GetReferenceData();
    const string header = bam.GetHeaderText();
    cout << "fix me!! headerSyasCoordSorted\n";
    return true;
}

}

/* namespace parser */
