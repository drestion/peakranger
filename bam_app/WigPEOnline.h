/*
 * WigPEOnline.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef WIGPEONLINE_H_
#define WIGPEONLINE_H_

#include "bam_app/OnlineBamMultiReportAppImp.h"
#include "utils/Tracer.h"
#include "bar/BarHist.h"
#include "bar/BarAux.h"
#include <algorithm>
namespace bam_app {
namespace aux {

class WigPEOnline: public OnlineBamMultiReportAppImp {
public:
    WigPEOnline(utils::TimeStampTracer& tracer);
    virtual ~WigPEOnline();
    void process(const BamTools::BamAlignment & read,
            const BamTools::RefVector & ref);
    void report(std::string& file_prefix);

    bool isGzip() const
    {
        return gz;
    }

    void setGzip(bool gzip)
    {
        this->gz = gzip;
    }

    uint32_t getNegCnt() const
    {
        return mNegCnt;
    }

    uint32_t getPosCnt() const
    {
        return mPosCnt;
    }

    bool isSplitByChr() const
    {
        return mbSplitByChr;
    }

    void setSplitByChr(bool mbSplitByChr)
    {
        this->mbSplitByChr = mbSplitByChr;
    }

    bool isSplitByStrand() const
    {
        return mbSplitByStrand;
    }

    void setSplitByStrand(bool mbSplitByStrand)
    {
        this->mbSplitByStrand = mbSplitByStrand;
    }

    uint32_t getExt() const
    {
        return mExt;
    }

    void setExt(uint32_t ext)
    {
        mExt = ext;
    }

private:
    std::map<std::string,ranger::bar::BarHist> mCounter;
    std::map<std::string,ranger::bar::BarHist> mNCounter;
    uint32_t mExt;
    bool mbSplitByStrand;
    bool mbSplitByChr;
    uint32_t mPosCnt;
    uint32_t mNegCnt;
    bool gz;
};

} /* namespace aux */
} /* namespace bam_app */
#endif /* WIGPEONLINE_H_ */
