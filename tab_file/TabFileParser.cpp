/*
 * TabFileParser.cpp
 *
 *  Created on: Mar 7, 2012
 *      Author: xinfeng
 */

#include "TabFileParser.h"
#include "utils/logger.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <assert.h>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;
using namespace boost;
namespace tab_file {

typedef size_t index_t;
const index_t TabFileParser::chr_col = 3 - 1;
const index_t TabFileParser::dir_col = 4 - 1;
const index_t TabFileParser::gene_start_col = 5 - 1;
const index_t TabFileParser::gene_end_col = 6 - 1;
const index_t TabFileParser::utr5_end_col = 7 - 1;
const index_t TabFileParser::utr3_start_col = 8 - 1;
const index_t TabFileParser::exon_start_col = 10 - 1;
const index_t TabFileParser::exon_end_col = 11 - 1;
const index_t TabFileParser::gene_name2 = 13 - 1;

TabFileParser::TabFileParser(const char* file,
        std::map<std::string, std::vector<TabGene> >& result) :
        _result(result), _file(file) {

}

TabFileParser::TabFileParser(std::string& file,
        std::map<std::string, std::vector<TabGene> >& result) :
        _result(result), _file(file) {

}

TabFileParser::~TabFileParser() {

}

void TabFileParser::parse() {
    using namespace boost;
    using namespace std;

    ifstream ifs(_file.c_str());

    string line;
    skipHeadLines(1, "#", ifs);
    while (getline(ifs, line)) {
        typedef tokenizer<char_separator<char> > tokenizer;
        char_separator<char> sep("\t");
        tokenizer tok(line, sep);

        vector<string> tokens;

        for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it) {
            tokens.push_back(*it);
        }

        parseTokens(tokens);
    }

    typedef std::map<std::string, std::vector<TabGene> > result_t;
    result_t::iterator it = _result.begin();
    for (; it != _result.end(); it++) {
        sort(it->second.begin(), it->second.end());
    }

}

void TabFileParser::parseTokens(std::vector<std::string>& tokens) {
    using namespace std;
    using namespace boost;
    typedef size_t loc_t;

    string chr;
    loc_t utr5_start, utr3_start, utr5_end, utr3_end;
    bool dir=true;

    vector<loc_t> exon_start, exon_end;

    chr = tokens[chr_col];

    if (tokens[dir_col] == "+") {
        dir = true;
    } else if (tokens[dir_col] == "-") {
        dir = false;
    } else {
        cout << "dir col is " << tokens[dir_col] << ". expecting + or -"
                << endl;
        throw std::runtime_error("bad_format");
    }

    utr5_start = lexical_cast < loc_t > (tokens[gene_start_col]);
    utr3_end = lexical_cast < loc_t > (tokens[gene_end_col]);
    utr5_end = lexical_cast < loc_t > (tokens[utr5_end_col]);
    utr3_start = lexical_cast < loc_t > (tokens[utr3_start_col]);
    tokenizer<> es(tokens[exon_start_col]);
    for (tokenizer<>::iterator it = es.begin(); it != es.end(); ++it) {
        exon_start.push_back(lexical_cast < loc_t > (*it));
    }
    tokenizer<> ee(tokens[exon_end_col]);
    for (tokenizer<>::iterator it = ee.begin(); it != ee.end(); ++it) {
        exon_end.push_back(lexical_cast < loc_t > (*it));
    }
    vector < Region > exons;
    for (size_t i = 0; i < exon_start.size(); i++) {
        exons.push_back(Region(exon_start[i], exon_end[i]));
    }
    assert(exons.size() > 0);
    exons.at(0).setL(utr5_end + 1);
    exons.at(exons.size() - 1).setR(utr3_start);
    if (utr5_start == utr5_end) {
        return;
    }
    if (utr5_end == utr3_start) {
        return;
    }
    TabGene gene(Region(utr5_start, utr5_end), Region(utr3_start, utr3_end),
            exons, tokens[gene_name2], dir);
    _result[chr].push_back(gene);
}



void TabFileParser::skipHeadLines(size_t n, const char* prefix,
        std::istream& is) {
    string line;
    while (n > 0 && getline(is, line)) {
        n--;
    }
}
}

