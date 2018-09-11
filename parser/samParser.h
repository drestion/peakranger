/*
 * samParser.h
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#ifndef SAMPARSER_H_
#define SAMPARSER_H_

#include "readsParser.h"
#include <iostream>
#include "short_reads/reads.h"


class samParser: public readsParser {
public:
    samParser() :
            readsParser() {

    }

    void parse(Reads& reads, std::string& filename);

    void parse(Reads& reads, std::istream& is);

    void parse(Reads& reads, std::istream& is,
            std::vector<std::string>& chrs_to_parse);

    void parse(Reads& reads, std::string& filename,
            std::vector<std::string>& chrs_to_parse);
};

#endif /* SAMPARSER_H_ */
