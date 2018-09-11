/*
 * bowtieParser.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#include "bowtieParser.h"
#include <string.h>
#include <istream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include <algorithm>
#include "utils/logger.h"
#include "utils/exceptions.h"
#include "short_reads/reads.h"
using namespace std;
//void bowtieParser::parse(Reads& outputreads) {
//	char _line[4096];
//	char line[4096];
//	char dir;
//	char cname[100];
//	char seq[4096];
//	int ret;
//	int cnum;
//	uint32_t start;
//	FILE* file = fopen(_file.c_str(),"r");
//	char BOWTIE_SCAN_FORMAT[496] = "%*[^\t]%*c%[^\t]%*c%[^\t]%*c%d%*c%[^\t]%*c";
//	while (fgets(_line, 4096, file) != NULL) {
//
//		strcpy(line, _line);
//		ret = sscanf(line, BOWTIE_SCAN_FORMAT, &dir, cname, &start, seq);
//
//		if (ret != 4) {
//			throw 4;
//		}
//
//		//Insert starts and increase window reads
//		if (dir == '-') {
//			start = start + 200;
//		}
//		string ga(cname);
//		outputreads.insertRead(ga, start);
//
//	}
//	fclose(file);
//}
//

/*
 * Performance data:
 * 10m reads
 scanf
 real	1m6.301s
 user	0m7.344s
 sys	0m31.170s

 istringstream
 real	1m10.789s
 user	0m15.133s
 sys	0m29.106s
 *
 */
void bowtieParser::parse(Reads& outputreads,
                         string& filename) {

    ifstream ifs(filename.c_str());

    if (!(ifs.good())) {
        throw FileNotGood(filename);
    }

    parse(outputreads,
          ifs);

    ifs.close();
}

void bowtieParser::parse(Reads& outputreads,
                         istream& ifs) {

    char dir;
    size_t ploc;
    string chr, seq, line;
    uint32_t loc;

    if (!(ifs.good())) {
        throw FileNotGood("input stream");
    }

    while (getline(ifs,
                   line)) {

        ploc = line.find("\t+\t");
        if (ploc == string::npos) {
            ploc = line.find("\t-\t");
            if (ploc == string::npos) {
                cout << "Warning: Didnt find direction character (+/-)"
                    " in the line: \n" << line << "\nThis line was skipped"
                << endl;
                continue;
            }
        }
        string nameTrimmed(line.begin() + ploc,
                           line.end());
        istringstream iss(nameTrimmed);
        iss >> dir >> chr >> loc >> seq;
        LOG_DEBUG4("In bowtie parser, Inserted read: "<<loc);
        if (dir == '-') {
            outputreads.neg_reads.insertRead(chr,
                                             loc);
        } else {
            outputreads.pos_reads.insertRead(chr,
                                             loc);
        }
    }
    outputreads.setReadlength(seq.size());
}

void bowtieParser::parse(Reads & outputreads,
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

void bowtieParser::parse(Reads & outputreads,
                         istream & ifs,
                         vector<string> & chrs_to_parse) {
    char dir;
    size_t ploc;
    string chr, seq, line;
    uint32_t loc;

    if (!(ifs.good())) {
        throw FileNotGood("The specified bowtie input file");
    }

    while (getline(ifs,
                   line)) {

        ploc = line.find("\t+\t");
        if (ploc == string::npos) {
            ploc = line.find("\t-\t");
            if (ploc == string::npos) {
                cout << "Warning: Didnt find direction character (+/-)"
                    " in the line: \n" << line << "\nThis line was skipped"
                << endl;
                continue;
            }
        }
        string nameTrimmed(line.begin() + ploc,
                           line.end());
        istringstream iss(nameTrimmed);
        iss >> dir >> chr >> loc >> seq;
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
            //Not the chr we want
        }
    }
    outputreads.setReadlength(seq.size());
}

