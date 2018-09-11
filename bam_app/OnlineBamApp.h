/*
 * OnlineBamApp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef ONLINEBAMAPP_H_
#define ONLINEBAMAPP_H_
#include <string>
#include <ostream>
#include <stdint.h>

#include "OnlineBamAppImp.h"
namespace bam_app {

/*
 * A framework for apps scale with
 * a single bam read using a single
 * bam file.
 */

class OnlineBamApp {
public:
    OnlineBamApp();
    virtual ~OnlineBamApp();

    OnlineBamApp(aux::OnlineBamAppImp* imp);

    virtual void processReads(const std::string& file, std::ostream& os);
    virtual void report(std::ostream& os);

    uint32_t getCntToReport() const
    {
        return mCntToReport;
    }

    void setCntToReport(uint32_t cntToReport)
    {
        mCntToReport = cntToReport;
    }

private:
    aux::OnlineBamAppImp* mImp;
//    uint32_t mCnt;
    uint32_t mCntToReport;
};

} /* namespace bam_app */
#endif /* ONLINEBAMAPP_H_ */
