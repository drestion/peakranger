/*
 * bedParser.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#include "bedParser.h"
#include "utils/logger.h"

#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include "utils/exceptions.h"
#include "short_reads/reads.h"
using namespace std;

void bedParser::parse(Reads & outputreads,
                      string & filename) {
    ifstream ifs(filename.c_str());

    if (!(ifs.good())) {
        throw FileNotGood(filename);
    }

    parse(outputreads,
          ifs);

    ifs.close();
}

void bedParser::parse(Reads & outputreads,
                      istream & ifs) {
    LOG_DEBUG1("Entering bedParser::parse");
    char dir;
    size_t ploc;
    string chr, seq, line;
    uint32_t loc, locend;

    if (!(ifs.good())) {
        throw FileNotGood("input stream");
    }

    while (getline(ifs,
                   line)) {

        ploc = line.find('+');
        if (ploc == string::npos) {
            ploc = line.find('-');
            if (ploc == string::npos) {
                cout << "Didnt find direction character (+/-)"
                " in the line: \n" << line
                << "\nThis line was skipped\n";
                continue;
            }
        }
        istringstream iss(line);
        iss >> chr >> loc >> locend >> seq >> seq >> dir;
        if(loc > locend) continue;
        if (dir == '-') {
            LOG_DEBUG5("INS: - "<<locend);
            outputreads.neg_reads.insertRead(chr,
                                             loc);//todo: CCAT requires locend
        } else {
            LOG_DEBUG5("INS: + "<<loc);
            outputreads.pos_reads.insertRead(chr,
                                             loc);
        }

    }

    outputreads.setReadlength(locend - loc);
    LOG_DEBUG1("Leaving bedParser::parse");
}

void bedParser::parse(Reads & outputreads,
                      string & filename,
                      vector<string> & chrs_to_parse) {
    ifstream ifs(filename.c_str());

    if (!(ifs.good())) {
        throw FileNotGood(filename);
    }

    parse(outputreads,
          ifs,
          chrs_to_parse);

    ifs.close();
}
void bedParser::parse(Reads & outputreads,
                      istream & ifs,
                      vector<string> & chrs_to_parse) {
    cout << "lala";
    char dir;
    size_t ploc;
    string chr, seq, line;
    uint32_t loc, locend;

    if (!(ifs.good())) {
        throw FileNotGood("The specified bowtie input file");
    }

    while (getline(ifs,
                   line)) {

        ploc = line.find('+');
        if (ploc == string::npos) {
            ploc = line.find('-');
            if (ploc == string::npos) {
                cout << "Didnt find direction character (+/-)"
                " in the line: \n" << line << "\nThis line was skipped" << endl;
                continue;
            }
        }
        istringstream iss(line);
        iss >> chr >> loc >> locend >> seq >> seq >> dir;

        if (std::find(chrs_to_parse.begin(),
                      chrs_to_parse.end(),
                      chr) != chrs_to_parse.end()) {
            LOG_DEBUG4("In bowtie parser, Inserted read: "<<loc);
            if (dir == '-') {
                outputreads.neg_reads.insertRead(chr,
                                                 loc);
            } else {
                outputreads.pos_reads.insertRead(chr,
                                                 loc);
            }
        } else {

        }
    }
    outputreads.setReadlength(seq.size());
}
