/*
 * NearbyGeneFinder.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: xfeng
 */

#include "tab_file/NearbyGeneFinder.h"
#include <algorithm>
using namespace std;
namespace tab_file {

NearbyGeneFinder::NearbyGeneFinder() :
        mSpan(0), mFile(""), mGenes(), mParser(), mParsed(false) {

}

NearbyGeneFinder::~NearbyGeneFinder() {

}

void print_gene_annot(const string& chr, ostream& os) {
}

void NearbyGeneFinder::removeDuplicates(vector<TabGene>& genes_to_filter) {

    sort(genes_to_filter.begin(), genes_to_filter.end(), name_lex_lt);

    genes_to_filter.erase(
            unique(genes_to_filter.begin(), genes_to_filter.end(), sameID),
            genes_to_filter.end());
}

void NearbyGeneFinder::setAnnoFile(const string& file) {
    mFile = file;
}

void NearbyGeneFinder::setSearchSpan(const uint32_t span) {
    mSpan = span;
}

void NearbyGeneFinder::getOverlappedGenes(const std::string& chr,
        const ranger::concepts::RegionUint32& region,
        vector<TabGene>& overlapedGenes) {
    if (!mParsed) {
        mParser.parse(mFile, mGenes);
        mParsed = true;
    }

    uint32_t left, right;
    left = region.getL() > mSpan ? region.getL() - mSpan : 0;
    right = region.getR() + mSpan;

    get_genes_to_filter(chr, left, right, mGenes, overlapedGenes);
    removeDuplicates(overlapedGenes);
}

void NearbyGeneFinder::get_genes_to_filter(const string& chr,
        const uint32_t left, const uint32_t right,
        std::map<std::string, std::vector<TabGene> >& parsedGenes,
        vector<TabGene>& genes_to_filter) {

    pair<vector<TabGene>::iterator, vector<TabGene>::iterator> bounds;
    TabGene tgg;
    tgg.utr5.setL(left);
    tgg.utr3.setR(right);
    bounds = equal_range(parsedGenes[chr].begin(), parsedGenes[chr].end(), tgg);
    copy(bounds.first, bounds.second, back_inserter(genes_to_filter));
}

} /* namespace tab_file */
