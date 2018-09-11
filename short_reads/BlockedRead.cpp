/*
 * BlockedRead.cpp
 *
 *  Created on: May 10, 2012
 *      Author: xin
 */

#include "BlockedRead.h"
#include <algorithm>

using namespace std;
namespace reads {
BlockedRead::BlockedRead() :
        mBlocks() {
}

BlockedRead::~BlockedRead() {
}

BlockedRead::BlockedRead(const vector<Read> & blocks) :
        mBlocks(blocks) {

}

void BlockedRead::getBlocks(vector<Read> & res) const {
    copy(mBlocks.begin(), mBlocks.end(), back_inserter(res));
}

std::vector<Read> BlockedRead::getBlocks() const {
    return mBlocks;
}

string BlockedRead::getChr() const {
    if (mBlocks.size())
        return mBlocks[0].getChr();
    else
        return "";
}

Strand BlockedRead::getDir() const {
    if (mBlocks.size()) {
        return mBlocks[0].getStrand();
    } else {
        return Strand();
    }
}

int32_t BlockedRead::getStart() const {
    if (mBlocks.size()) {
        return mBlocks[0].getStart();
    } else {
        return 0;
    }
}

}

void reads::offset(BlockedRead& read, int32_t length) {
    vector<Read> rds;
    read.getBlocks(rds);
    foreach(Read& r, rds) {
        r.setStart(r.getStart() + length);
        r.setEnd(r.getEnd() + length);
    }
    read.setBlocks(rds);
}

