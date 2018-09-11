/*
 * UCSCAnnoParser.h
 *
 *  Created on: Jun 29, 2012
 *      Author: xfeng
 */

#ifndef UCSCANNOPARSER_H_
#define UCSCANNOPARSER_H_
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <stdint.h>
#include "TabGene.h"
namespace tab_file {

class UCSCAnnoParser {
public:
    UCSCAnnoParser();
    virtual ~UCSCAnnoParser();
    void parse(std::string& file,
            std::map<std::string, std::vector<TabGene> >& result);

private:
    void parseTokens(std::vector<std::string>& tokens,
            std::map<std::string, std::vector<TabGene> >& result);
    void skipHeadLines(size_t n, const char* prefix, std::istream& is);

private:

    typedef size_t index_t;

    static const index_t chr_col;
    static const index_t dir_col;
    static const index_t gene_start_col;
    static const index_t gene_end_col;
    static const index_t utr5_end_col;
    static const index_t utr3_start_col;
    static const index_t exon_start_col;
    static const index_t exon_end_col;
    static const index_t gene_name2;
};

} /* namespace tab_file */
#endif /* UCSCANNOPARSER_H_ */
