/*
 * BamBlockPEReadsParserhotelImp.h
 *
 *  Created on: May 13, 2012
 *      Author: xfeng
 */

#ifndef BAMBLOCKPEREADSPARSERHOTELIMP_H_
#define BAMBLOCKPEREADSPARSERHOTELIMP_H_

#include "bamtools/BamAux.h"
#include "short_reads/BlockedRead.h"
#include "short_reads/PairEndedReads.h"
#include "BamBlockPEReadsParserImp.h"
namespace parser {
namespace aux {

class BamBlockPEReadsParser_hotelImp:public BamBlockPEReadsParserImp {
public:
    BamBlockPEReadsParser_hotelImp();
    virtual ~BamBlockPEReadsParser_hotelImp();

    void parse(const BamTools::BamAlignment & read,
            reads::PairEndedReads<reads::BlockedRead> & reads,
            const BamTools::RefVector & ref);
    void flush(
               reads::PairEndedReads<reads::BlockedRead> & reads,const BamTools::RefVector& ref);


    bool checkedIn(const BamTools::BamAlignment& read);

    void checkIn(const BamTools::BamAlignment& bam);
    BamTools::BamAlignment checkOut(const BamTools::BamAlignment& bam);


private:
    std::vector<BamTools::BamAlignment> mUnMatched;
};

} /* namespace aux */
} /* namespace parser */
#endif /* BAMBLOCKPEREADSPARSERHOTELIMP_H_ */
