/*
 * BamHit.h
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#ifndef BAMHIT_H_
#define BAMHIT_H_

#include "bamtools/BamAux.h"
#include "Read.h"
namespace reads {

class BamHit: public Read {
public:
    BamHit();
    BamHit(const BamTools::BamAlignment& bam, const BamTools::RefVector& ref);
    BamHit(const BamTools::BamAlignment& bam, const std::string& chr);
    BamHit(const BamTools::BamAlignment& bam, const char* chr);
    virtual ~BamHit();


private:

    void getStartEndDir(const BamTools::BamAlignment& bam);
    int32_t getReadLength(const BamTools::BamAlignment& bam);

};

} /* namespace reads */
#endif /* BAMHIT_H_ */
