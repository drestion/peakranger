/*
 * SimpleThresholder.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef SIMPLETHRESHOLDER_H_
#define SIMPLETHRESHOLDER_H_

#include "StockBamMultipleDatasetsAppImp.h"
#include "bar/RegionCount.h"
#include "utils/Tracer.h"
namespace bam_app {
namespace aux {

class SimpleThresholder: public bam_app::aux::StockBamMultipleDatasetsAppImp {
public:
    SimpleThresholder(utils::TimeStampTracer& tracer);
    virtual ~SimpleThresholder();
    void report(
             std::vector<reads::PairEndedReads<reads::BlockedRead> >& rdsVec,
            std::ostream& os);


    int32_t getLowerCutoff() const;
    void setLowerCutoff(int32_t lowerCutoff);
    int32_t getUpperCufoff() const;
    void setUpperCufoff(int32_t upperCufoff);

private:
    int32_t mLowerCutoff;
    int32_t mUpperCufoff;
    bool isInCutRegion(ranger::bar::RegionCount & r, int32_t lhs, int32_t rhs);
};

} /* namespace aux */
} /* namespace bam_app */
#endif /* SIMPLETHRESHOLDER_H_ */
