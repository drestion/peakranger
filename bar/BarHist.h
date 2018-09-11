/*
 * BarHist.h
 *
 *  Created on: May 24, 2012
 *      Author: tania
 */

#ifndef BARHIST_H_
#define BARHIST_H_

#include "BarCounter.h"
#include "concepts/RegionInt32.h"
#include "RegionCount.h"
#include "short_reads/Read.h"
#include "bamtools/BamAux.h"
#include <vector>
namespace ranger {
namespace bar {

class BarHist {
    friend void operator+=(BarHist& b1, const BarHist& b2);
    friend void operator+=(std::map<std::string, BarHist>& b1,
             std::map<std::string, BarHist>& b2);
    reads::Read extendReadToLength(const reads::Read& read, int32_t length);
public:
    BarHist();
    virtual ~BarHist();

    void add(const reads::Read& read);
    void add(const BamTools::BamAlignment& bam);
    void add(const std::vector<reads::Read>& hits);
    void getCount(const std::vector<reads::Read>& hits,
            std::vector<RegionCount>& res);
    void reset();

    const boost::icl::BarCounter& counter() const {
        return mCounter;
    }

private:
    boost::icl::BarCounter mCounter;

};

} /* namespace bar */
}

/* namespace ranger */
#endif /* BARHIST_H_ */
