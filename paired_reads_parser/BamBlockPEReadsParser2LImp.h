/*
 * BamBlockPEReadsParser2LImp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef BAMBLOCKPEREADSPARSER2LIMP_H_
#define BAMBLOCKPEREADSPARSER2LIMP_H_

#include "BamBlockPEReadsParserImp.h"

namespace parser {
namespace aux {

/*
 * This parser works as long as the bam file is sorted
 * by read name or chromosome.
 * The parser assumes two adjacent reads have the same
 * read name so they form a read pair. This is true when
 * the bam file is sorted by read name. (-n in samtools sort)
 * When this assumption does not hold, the parser will
 * still be able to make pseudo pairs out of two adjacent reads.
 * It depends on the user to decide if this behavior is
 * OK.
 */

class BamBlockPEReadsParser2LImp: public parser::aux::BamBlockPEReadsParserImp {
public:
    BamBlockPEReadsParser2LImp();
    virtual ~BamBlockPEReadsParser2LImp();
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
#endif /* BAMBLOCKPEREADSPARSER2LIMP_H_ */
