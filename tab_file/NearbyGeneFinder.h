/*
 * NearbyGeneFinder.h
 *
 *  Created on: Jun 29, 2012
 *      Author: xfeng
 */

#ifndef NEARBYGENEFINDER_H_
#define NEARBYGENEFINDER_H_

#include <vector>
#include <stdint.h>
#include <string>
#include "tab_file/TabGene.h"
#include "concepts/RegionInt32.h"
#include "tab_file/UCSCAnnoParser.h"
namespace tab_file {

class NearbyGeneFinder {
public:
    NearbyGeneFinder();
    virtual ~NearbyGeneFinder();

    void setAnnoFile(const std::string& file);
    void setSearchSpan(const uint32_t span);
    void getOverlappedGenes(const std::string& chr,
            const ranger::concepts::RegionUint32& region,
            std::vector<TabGene>& overlaped);

private:
    void removeDuplicates(std::vector<tab_file::TabGene>& genes_to_filter);
    void get_genes_to_filter(const std::string& chr, const uint32_t left,
            const uint32_t right,
            std::map<std::string, std::vector<TabGene> >& parsedGenes,
            std::vector<tab_file::TabGene>& genes_to_filter);

private:

    uint32_t mSpan;
    std::string mFile;
    std::map<std::string, std::vector<TabGene> > mGenes;
    UCSCAnnoParser mParser;
    bool mParsed;
};

} /* namespace tab_file */
#endif /* NEARBYGENEFINDER_H_ */
