/*
 * chrt.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef CHRT_H_
#define CHRT_H_

#include "common/stl_header.h"
#include "ccat_peak.h"
namespace ccat_aux {

class chr_t {
    friend std::ostream& operator<<(std::ostream& os, const chr_t& chr) {
        os << chr.chromName << "\tl1PosTAGS[0] " << chr.l1PosTags[0]
                << "\tl1PosTAGS[last] "
                << chr.l1PosTags[chr.l1PosTags.size() - 1] << "\tchr Size "
                << chr.chromSize << "\tl1NegTAGS[0] " << chr.l1NegTags[0]
                << "\tl1NegTAGS[last] "
                << chr.l1NegTags[chr.l1NegTags.size() - 1]
                << "\nl1PosTags.size() " << chr.l1PosTags.size()
                << "\tl1NegTags.size() " << chr.l1NegTags.size()
                << "\tl2PosTags.size() " << chr.l2PosTags.size()
                << "\tl2NegTags.size() " << chr.l2NegTags.size();
        return os;
    }
public:
    chr_t();
    virtual ~chr_t();

    std::string chromName;
    int chromIndex;
    size_t chromSize;

    std::vector<size_t> l1PosTags;
    std::vector<size_t> l1NegTags;
    std::vector<size_t> l2PosTags;
    std::vector<size_t> l2NegTags;
    std::vector<peak_t> l1Peaks;
    std::vector<peak_t> l2Peaks;

};

} /* namespace ccat_aux */
#endif /* CHRT_H_ */
