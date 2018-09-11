/*
 * bowtieParser.h
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#ifndef BOWTIEPARSER_H_
#define BOWTIEPARSER_H_

#include "readsParser.h"

class bowtieParser: public readsParser {
public:
    bowtieParser() :
            readsParser() {

    }
    void parse(Reads& reads, std::string& filename);

    void parse(Reads& reads, std::istream& is);

    void parse(Reads& reads, std::istream& is,
            std::vector<std::string>& chrs_to_parse);

    void parse(Reads& reads, std::string& filename,
            std::vector<std::string>& chrs_to_parse);

protected:
    void inline _parse_line(Reads& reads, std::string& line);
};

#endif /* BOWTIEPARSER_H_ */
