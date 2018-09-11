/*
 * BamBlockPEReadsParserImp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef BAMBLOCKPEREADSPARSERIMP_H_
#define BAMBLOCKPEREADSPARSERIMP_H_
#include "bamtools/BamAux.h"
#include "short_reads/BlockedRead.h"
#include "short_reads/PairEndedReads.h"

namespace parser {
namespace aux {

class BamBlockPEReadsParserImp {
public:
    BamBlockPEReadsParserImp();
    virtual ~BamBlockPEReadsParserImp();

    virtual void parse(const BamTools::BamAlignment & read,
            reads::PairEndedReads<reads::BlockedRead> & reads,
            const BamTools::RefVector & ref) ;

    virtual void flush(reads::PairEndedReads<reads::BlockedRead> & reads,
            const BamTools::RefVector& ref);

    std::string getParserId() const
    {
        return parserID;
    }

    void setParserId(std::string parserId)
    {
        parserID = parserId;
    }

protected:
    std::string parserID;
};

} /* namespace aux */
} /* namespace parser */
#endif /* BAMBLOCKPEREADSPARSERIMP_H_ */
