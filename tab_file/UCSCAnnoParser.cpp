/*
 * UCSCAnnoParser.cpp
 *
 *  Created on: Jun 29, 2012
 *      Author: xfeng
 */

#include "tab_file/UCSCAnnoParser.h"

#include "utils/logger.h"
#include "utils/Guarded.h"
#include "utils/exceptions.h"
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
const index_t UCSCAnnoParser::chr_col = 3 - 1;
const index_t UCSCAnnoParser::dir_col = 4 - 1;
const index_t UCSCAnnoParser::gene_start_col = 5 - 1;
const index_t UCSCAnnoParser::gene_end_col = 6 - 1;
const index_t UCSCAnnoParser::utr5_end_col = 7 - 1;
const index_t UCSCAnnoParser::utr3_start_col = 8 - 1;
const index_t UCSCAnnoParser::exon_start_col = 10 - 1;
const index_t UCSCAnnoParser::exon_end_col = 11 - 1;
const index_t UCSCAnnoParser::gene_name2 = 13 - 1;
UCSCAnnoParser::UCSCAnnoParser() {

}

UCSCAnnoParser::~UCSCAnnoParser() {
}

void UCSCAnnoParser::parse(std::string& f,
        std::map<std::string, std::vector<TabGene> >& result) {

    ifstream ifs(f.c_str());
    utils::Guarded<FileNotGood> g(!ifs, f.c_str());
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
        if (gene_name2 > tokens.size()) {
            //assuming gene_name2 is the larget one!
            continue;
        }
        parseTokens(tokens, result);
    }

    typedef std::map<std::string, std::vector<TabGene> > result_t;
    result_t::iterator it = result.begin();
    for (; it != result.end(); it++) {
        sort(it->second.begin(), it->second.end());
    }
}

void UCSCAnnoParser::parseTokens(std::vector<std::string>& tokens,
        std::map<std::string, std::vector<TabGene> >& _result) {

    typedef size_t loc_t;

    string chr;
    loc_t utr5_start, utr3_start, utr5_end, utr3_end;
    bool dir = true;

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

    utr5_start = lexical_cast<loc_t>(tokens[gene_start_col]);
    utr3_end = lexical_cast<loc_t>(tokens[gene_end_col]);
    utr5_end = lexical_cast<loc_t>(tokens[utr5_end_col]);
    utr3_start = lexical_cast<loc_t>(tokens[utr3_start_col]);
    tokenizer<> es(tokens[exon_start_col]);
    for (tokenizer<>::iterator it = es.begin(); it != es.end(); ++it) {
        exon_start.push_back(lexical_cast<loc_t>(*it));
    }
    tokenizer<> ee(tokens[exon_end_col]);
    for (tokenizer<>::iterator it = ee.begin(); it != ee.end(); ++it) {
        exon_end.push_back(lexical_cast<loc_t>(*it));
    }
    vector<Region> exons;
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

void UCSCAnnoParser::skipHeadLines(size_t n, const char* prefix,
        std::istream& is) {
    string line;
    while (n > 0 && getline(is, line)) {
        n--;
    }
}
}

