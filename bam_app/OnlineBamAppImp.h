/*
 * OnlineBamAppImp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef ONLINEBAMAPPIMP_H_
#define ONLINEBAMAPPIMP_H_
#include "bamtools/BamAux.h"
#include <ostream>
namespace bam_app {
namespace aux {

class OnlineBamAppImp {
public:
    OnlineBamAppImp();
    virtual ~OnlineBamAppImp();
    virtual void process(const BamTools::BamAlignment & read,
            const BamTools::RefVector & ref);
    virtual void report(std::ostream& os);
};

} /* namespace aux */
} /* namespace bam_app */
#endif /* ONLINEBAMAPPIMP_H_ */
