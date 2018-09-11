/*
 * bamParser.h
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#ifndef BAMPARSER_H_
#define BAMPARSER_H_

#include "readsParser.h"
#include "bamtools/BamAux.h"

class bamParser: public readsParser {

public:
    void parse(Reads& reads, std::string& filename);

    void parse(Reads& reads, std::istream& is);

    void parse(Reads& reads, std::istream& is,
            std::vector<std::string>& chrs_to_parse);

    void parse(Reads& reads, std::string& filename,
            std::vector<std::string>& chrs_to_parse);

private:

    void insertRead(const BamTools::BamAlignment& read, Reads& reads,
            std::string& chr);
    void updateAvgReadLength(uint64_t& readCnt, uint32_t& meanReadLen,
            BamTools::BamAlignment& read);

};

#endif /* BAMPARSER_H_ */
