/*
 * readsParser.h
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#ifndef READSPARSER_H_
#define READSPARSER_H_

#include <string>
#include <iostream>
#include <vector>
class Reads;
class readsParser {
public:
    readsParser() {

    }

    virtual ~readsParser() {
    }

public:
    void virtual parse(Reads& reads, std::string& filename);

    void virtual parse(Reads& reads, std::istream& is);

    void virtual parse(Reads& reads, std::istream& is,
            std::vector<std::string>& chrs_to_parse);

    void virtual parse(Reads& reads, std::string& filename,
            std::vector<std::string>& chrs_to_parse);
};

#endif /* READSPARSER_H_ */
