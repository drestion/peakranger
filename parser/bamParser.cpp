/*
 * bamParser.cpp
 *
 *  Created on: Apr 29, 2011
 *      Author: xin
 */

#include "bamParser.h"
#include "short_reads/bamFormatAux.h"
#include "bamtools/BamReader.h"
#include "utils/exceptions.h"
#include "common/ranger_debug.h"
#include <string.h>
#include <istream>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdint.h>
#include <algorithm>

using namespace std;
using namespace reads;

void bamParser::insertRead(const BamTools::BamAlignment& read, Reads& reads,
		string& chr) {
	int32_t loc = read.Position;
	bool dir;

	dir = (read.IsReverseStrand() ? false : true);
	if (loc > 0) {
		uint32_t tmp = (uint32_t) loc;
		if (dir) {
			reads.pos_reads.insertRead(chr, tmp);
		} else {
			reads.neg_reads.insertRead(chr, tmp);
		}
	}
}

void bamParser::updateAvgReadLength(uint64_t& readCnt, uint32_t& meanReadLen,
		BamTools::BamAlignment& read) {
	readCnt++;


	uint64_t currentMeanLen = ( reads::getReadLength(read)
			+ meanReadLen * (readCnt - 1)) / readCnt;
	meanReadLen = (uint32_t) currentMeanLen;
}

void bamParser::parse(Reads & reads, string & filename,
		vector<string> & chrs_to_parse) {
	BamTools::BamReader bam;
	BamTools::BamAlignment read;

	string chr;
	uint64_t readCnt = 0;
	uint32_t meanReadLen = 0;

	if (!(bam.Open(filename))) {
		throw FileNotGood(filename);
	}

	const BamTools::RefVector refvec = bam.GetReferenceData();

	while (bam.GetNextAlignment(read)) {
		chr = getR1Chr(read, refvec);
		if (isGoodRead(read)) {
			if (isChrToParse(chrs_to_parse, chr)) {
				updateAvgReadLength(readCnt, meanReadLen, read);
				insertRead(read, reads, chr);
			}
		}
	}
	reads.setReadlength(meanReadLen);
}

void bamParser::parse(Reads & reads, string & filename) {
	BamTools::BamReader bam;
	BamTools::BamAlignment read, pre_read;
	string chr;
	uint64_t readCnt = 0;
	uint32_t meanReadLen = 0;

	if (!(bam.Open(filename))) {
		throw FileNotGood(filename);
	}
	const BamTools::RefVector refvec = bam.GetReferenceData();
	while (bam.GetNextAlignment(read)) {
		chr = getR1Chr(read, refvec);
		if (isGoodRead(read)) {
			updateAvgReadLength(readCnt, meanReadLen, read);
			insertRead(read, reads, chr);
		}
	}
	reads.setReadlength(meanReadLen);
}

void bamParser::parse(Reads & reads, istream & is,
		vector<string> & chrs_to_parse) {
	throw RangerException(
			"This bam parser has not been fully implemented yet.");
}

void bamParser::parse(Reads & reads, istream & is) {

	throw RangerException(
			"This bam parser has not been fully implemented yet.");
}
