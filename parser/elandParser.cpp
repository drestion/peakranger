/*
 * elandParser.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#include "elandParser.h"
//
//elandParser::elandParser() {
//	// TODO Auto-generated constructor stub
//
//}
//
//elandParser::~elandParser() {
//	// TODO Auto-generated destructor stub
//}
//
//void elandParser::parse(Reads& outputreads) {
//	char _line[FILENAME_LENGTH];
//	char line[FILENAME_LENGTH];
//	char dir;
//	char type[CHR_NAME_LENGTH];
//	char cname[CHR_NAME_LENGTH];
//	char seq[4096];
//	int ret;
//	int cnum;
//	map<string, int>::iterator it;
//	Dataline _result;
//
//	while (fgets(_line, FILENAME_LENGTH, file) != NULL) {
//
//		strcpy(line, _line);
//		ret = sscanf(line, ELAND_SCAN_FORMAT, seq, type, cname, &_result.start, &dir);
//
//		if (ret != 5) {
//			throw BadTextLineException(line);
//		}
//
//		if (!starts_with(type, "U", 1)) {
//			if (!starts_with(type, "u", 1)) {
//				continue; //Only accepts unique match.
//			}
//		}
//
//		_result.chr = cname;
//		if (dir == 'F') {
//			_result.orientation = true;
//		} else {
//			_result.orientation = false;
//		}
//		_result.seq = seq;
//
//		it = chr_ind.find(_result.chr);
//
//		if (it == chr_ind.end())
//			throw BadChrException(_result.chr);
//
//		cnum = getCnum(const_cast<char *> (_result.chr.c_str()));
//		if (!_result.orientation) {
//			start = _result.start + (_result.seq.size() - READ_LENGTH);
//
//		} else {
//			start = _result.start;
//
//		}
//		outputreads.insertRead(chr, start);
//	}
//}
