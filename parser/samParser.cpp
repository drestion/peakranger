/*
 * samParser.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#include "samParser.h"
#include "utils/exceptions.h"
#include "utils/logger.h"
#include "parser/sam_read.h"
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <utility>
using namespace std;

void samParser::parse(Reads& outputreads,
                      string& filename) {
    ifstream ifs(filename.c_str());
    if (!(ifs.good())) {
        throw FileNotGood(filename);
    }
    parse(outputreads,
          ifs);
    ifs.close();
}

void samParser::parse(Reads& outputreads,
                      istream& ifs) {
    bool dir, dir2;
    size_t ploc;
    string chr, seq, line, preLine, chr2, seq2;
    int32_t loc, loc2;

    map<string, string> unmapped;
    map<string, string>::iterator itr;
    sam_read read;
    while ((getline(ifs,
                    line))) {

        ploc = line.find('\t');
        if (ploc == string::npos) {
            cout << "Warning: invalid line skipped:\n" << line << endl;
            continue;
        }
        try {
            read.parseLine(line);
        }
        catch (DataLineNotValid& e) {
            cout << "Warning: can not parse the line:\n" << line << endl;
            continue;
        }
        chr = read.RefName;
        if (chr == "*") continue;
        if (chr == "") continue;
        if (chr == " ") continue;
        loc = read.Position;
        if (loc < 0) continue;
        seq = read.QueryBases;
        uint32_t tmp = (uint32_t) loc;
        dir = (read.IsReverseStrand() ? false : true);
        if (read.IsPaired() && read.IsProperPair() && !(read.IsFailedQC())
        && read.IsMateMapped() && read.IsFirstMate() && !(read.IsDuplicate())) {
            //paired reads
            loc2 = read.MatePosition;
            if (loc2 < 0) continue;
            dir2 = (read.IsMateReverseStrand() ? false : true);
            //keep only the start location
            if (dir && dir2) {
                loc = loc2 > loc ? loc : loc2;
                tmp = (uint32_t) loc;
                outputreads.pos_reads.insertRead(chr,
                                                 tmp);
            } else if (!dir && !dir2) {
                loc = loc2 > loc ? loc : loc2;
                tmp = (uint32_t) loc;
                outputreads.neg_reads.insertRead(chr,
                                                 tmp);
            } else {
                //treat the pair as two reads
                if (dir) {
                    tmp = (uint32_t) loc;
                    outputreads.pos_reads.insertRead(chr,
                                                     tmp);
                } else {
                    tmp = (uint32_t) loc;
                    outputreads.neg_reads.insertRead(chr,
                                                     tmp);
                }
                if (dir2) {
                    tmp = (uint32_t) loc2;
                    outputreads.pos_reads.insertRead(chr,
                                                     tmp);
                } else {
                    tmp = (uint32_t) loc2;
                    outputreads.neg_reads.insertRead(chr,
                                                     tmp);
                }
            }

        } else if (!(read.IsPaired()) && !(read.IsFailedQC())) {
            //unpaired reads
            if (dir) {
                outputreads.pos_reads.insertRead(chr,
                                                 tmp);
            } else {
                outputreads.neg_reads.insertRead(chr,
                                                 tmp);
            }
        } else {
            //skip all other reads
            // such as second mates.
        }
    }
    //assumes the same read length
    outputreads.setReadlength(seq.size());
}

void samParser::parse(Reads & outputreads,
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

void samParser::parse(Reads & outputreads,
                      istream & ifs,
                      vector<string> & chrs_to_parse) {
    bool dir, dir2;
    size_t ploc;
    string chr, seq, line, preLine, chr2, seq2;
    int32_t loc, loc2;

    map<string, string> unmapped;
    map<string, string>::iterator itr;
    sam_read read;
    while ((getline(ifs,
                    line))) {

        ploc = line.find('\t');
        if (ploc == string::npos) {
            cout << "Warning: invalid line skipped:\n" << line << endl;
            continue;
        }
        try {
            read.parseLine(line);
        }
        catch (DataLineNotValid& e) {
            cout << "Warning: invalid line skipped:\n" << line << endl;
            continue;
        }
        chr = read.RefName;
        if (chr == "*") continue;
        if (chr == "") continue;
        if (chr == " ") continue;
        loc = read.Position;
        if (loc < 0) continue;
        seq = read.QueryBases;
        uint32_t tmp = (uint32_t) loc;
        dir = (read.IsReverseStrand() ? false : true);
        if (std::find(chrs_to_parse.begin(),
                      chrs_to_parse.end(),
                      chr) != chrs_to_parse.end()) {
            if (read.IsPaired() && read.IsProperPair() && !(read.IsFailedQC())
            && read.IsMateMapped() && read.IsFirstMate()
            && !(read.IsDuplicate())) {
                //paired reads
                loc2 = read.MatePosition;
                if (loc2 < 0) continue;
                dir2 = (read.IsMateReverseStrand() ? false : true);
                //keep only the start location
                if (dir && dir2) {
                    loc = loc2 > loc ? loc : loc2;
                    tmp = (uint32_t) loc;
                    outputreads.pos_reads.insertRead(chr,
                                                     tmp);
                } else if (!dir && !dir2) {
                    loc = loc2 > loc ? loc : loc2;
                    tmp = (uint32_t) loc;
                    outputreads.neg_reads.insertRead(chr,
                                                     tmp);
                } else {
                    //treat the pair as two reads
                    if (dir) {
                        tmp = (uint32_t) loc;
                        outputreads.pos_reads.insertRead(chr,
                                                         tmp);
                    } else {
                        tmp = (uint32_t) loc;
                        outputreads.neg_reads.insertRead(chr,
                                                         tmp);
                    }
                    if (dir2) {
                        tmp = (uint32_t) loc2;
                        outputreads.pos_reads.insertRead(chr,
                                                         tmp);
                    } else {
                        tmp = (uint32_t) loc2;
                        outputreads.neg_reads.insertRead(chr,
                                                         tmp);
                    }
                }

            } else if (!(read.IsPaired()) && !(read.IsFailedQC())) {
                //unpaired reads
                if (dir) {
                    outputreads.pos_reads.insertRead(chr,
                                                     tmp);
                } else {
                    outputreads.neg_reads.insertRead(chr,
                                                     tmp);
                }
            } else {
                //skip all other reads
                // such as second mates.
            }
        }
    }
    //assumes the same read length
    outputreads.setReadlength(seq.size());
}
