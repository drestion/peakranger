/*
 * bedParser.h
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#ifndef BEDPARSER_H_
#define BEDPARSER_H_

#include "readsParser.h"

class bedParser: public readsParser {
public:
    bedParser() :
            readsParser() {
    }
    ~bedParser() {
    }
    void parse(Reads& reads, std::string& filename);

    void parse(Reads& reads, std::istream& is);

    void parse(Reads& reads, std::istream& is,
            std::vector<std::string>& chrs_to_parse);

    void parse(Reads& reads, std::string& filename,
            std::vector<std::string>& chrs_to_parse);

};

#endif /* BEDPARSER_H_ */
