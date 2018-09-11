/*
 * LibraryComplexity.h
 *
 *  Created on: Jul 10, 2012
 *      Author: xfeng
 */

#ifndef LIBRARYCOMPLEXITY_H_
#define LIBRARYCOMPLEXITY_H_

#include "StockBamMultipleDatasetsAppImp.h"
#include "utils/Tracer.h"
namespace bam_app {
namespace aux {

class LibraryComplexity: public StockBamMultipleDatasetsAppImp {
public:
    LibraryComplexity(utils::TimeStampTracer& tra);
    virtual ~LibraryComplexity();

    void report(std::vector<reads::PairEndedReads<reads::BlockedRead> >& rdsVec,
            std::ostream& os);
    bool mIsPE;
};

} /* namespace aux */
} /* namespace bam_app */
#endif /* LIBRARYCOMPLEXITY_H_ */
