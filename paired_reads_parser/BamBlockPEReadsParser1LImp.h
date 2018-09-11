/*
 * BamBlockPEReadsParser1LImp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef BAMBLOCKPEREADSPARSER1LIMP_H_
#define BAMBLOCKPEREADSPARSER1LIMP_H_

#include "BamBlockPEReadsParserImp.h"

namespace parser {
namespace aux {

class BamBlockPEReadsParser1LImp: public parser::aux::BamBlockPEReadsParserImp {
public:
    BamBlockPEReadsParser1LImp();
    virtual ~BamBlockPEReadsParser1LImp();
    void parse(const BamTools::BamAlignment & read,
            reads::PairEndedReads<reads::BlockedRead> & reads,
            const BamTools::RefVector & ref);

    void flush(reads::PairEndedReads<reads::BlockedRead> & reads,
            const BamTools::RefVector& ref);


private:
    BamTools::BamAlignment mPRead;
};

} /* namespace aux */
} /* namespace parser */
#endif /* BAMBLOCKPEREADSPARSER1LIMP_H_ */
