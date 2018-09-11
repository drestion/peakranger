/*
 * Reads.cpp
 *
 *  Created on: May 5, 2011
 *      Author: xin
 */

#include "reads.h"

#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include <algorithm>

#define MAX_CHR_NUM 10000
using namespace std;
typedef map<string, vector<uint32_t> > reads_t;
void Reads::_insertRead(string& chr, uint32_t read, reads_t& reads) {
    reads_t::iterator ind;
    rt_assert_neq_msg(chr, "", "The name of chr is empty.")
    rt_assert_neq_msg(chr, "*", "The name of chr is *.")
    if ((ind = reads.find(chr)) == reads.end()) {
        vector<uint32_t> tmp;
        reads[chr] = tmp;
    }
    reads[chr].push_back(read);
}
/*
 * return the number of pos and neg reads
 */
uint64_t Reads::size() const {
    return pos_reads.size() + neg_reads.size();
}

uint32_t Reads::getReadlength() const {
    return _readlength;
}

void Reads::setReadlength(uint32_t _readlength) {
    this->_readlength = _readlength;
}

pritr Reads::pos_reads::begin() const {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();
    return _reads._posReads.begin();
}

pritr Reads::pos_reads::end() const {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();
    return _reads._posReads.end();
}

ritr Reads::pos_reads::end_of(string & chr) const {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();
    if (hasReadsOnChr(chr)) {
        return _reads._posReads[chr].end();
    } else {
        string str("Chromosome ");
        str += chr;
        str += " was not found in parsed positive reads.";
        throw ChrNotFound(str);
    }
}

ritr Reads::pos_reads::begin_of(string & chr) const {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();
    if (hasReadsOnChr(chr)) {
        return _reads._posReads[chr].begin();
    } else {
        string str("Chromosome ");
        str += chr;
        str += " was not found in parsed positive reads.";
        throw ChrNotFound(str);
    }
}

void Reads::pos_reads::insertRead(string & chr, uint32_t & read) {
    _reads._insertRead(chr, read, _reads._posReads);
}

uint64_t Reads::pos_reads::size() const {
    uint64_t size = 0;
    pair<string, vector<uint32_t> > reads;
    foreach(reads,
            _reads._posReads) {
        size += (uint64_t) (reads.second.size());
    }
    return size;
}

vector<string> Reads::pos_reads::chrs() const {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();
    return _reads._posChrs;
}

bool Reads::pos_reads::hasReadsOnChr(string & chr) const {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();

    if (binary_search(_reads._posChrs.begin(), _reads._posChrs.end(), chr)) {
        return true;
    }
    return false;
}

void Reads::pos_reads::remove(string & chr) {
    if (!_reads._noMorePosReads)
        _reads.pos_reads._posReadsInsertionComplete();
    if (hasReadsOnChr(chr)) {
        _reads._posReads.erase(chr);
        vector<string>::iterator it;
        it = find(_reads._posChrs.begin(), _reads._posChrs.end(), chr);
        _reads._posChrs.erase(it);
    }
}

void Reads::pos_reads::_posReadsInsertionComplete() {
    _reads._noMorePosReads = true;
    _reads.extractChrs(_reads._posReads, _reads._posChrs);
    LOG_DEBUG1("Pos reads insertion complete.");
    _reads.sortReadsOnEachChromosome(_reads._posReads);
    _reads._pos_is_sorted = true;

}
/*
 * Access to reads in a region
 */
void Reads::pos_reads::getReads(string& chr, uint32_t start, uint32_t end,
        vector<uint32_t>::iterator& rstart, vector<uint32_t>::iterator& rend) {
    assert_gt(end, start)
    assert_neq(chr, "")

    vector<uint32_t>::iterator preadsstart, preadsend;
    if (hasReadsOnChr(chr)) {
        preadsstart = begin_of(chr);
        preadsend = end_of(chr);
        rstart = lower_bound(preadsstart, preadsend, start);
        rend = upper_bound(rstart, preadsend, end);
    } else {
        string str("Chromosome ");
        str += chr;
        str += " contains no positive reads.";
        throw ChrNotFound(str);
    }
}
/*
 * For negative reads
 */
pritr Reads::neg_reads::begin() const {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();
    return _reads._negReads.begin();
}

vector<string> Reads::neg_reads::chrs() const {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();
    return _reads._negChrs;
}

bool Reads::neg_reads::hasReadsOnChr(string & chr) const {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();

    if (binary_search(_reads._negChrs.begin(), _reads._negChrs.end(), chr)) {
        return true;
    }
    return false;
}
/*
 * Remove reads on chr if it has reads on it.
 */
void Reads::neg_reads::remove(string & chr) {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();
    if (hasReadsOnChr(chr)) {
        _reads._negReads.erase(chr);
        vector<string>::iterator it;
        it = find(_reads._negChrs.begin(), _reads._negChrs.end(), chr);
        _reads._negChrs.erase(it);
    }
}

ritr Reads::neg_reads::end_of(string & chr) const {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();
    if (hasReadsOnChr(chr)) {
        return _reads._negReads[chr].end();
    } else {
        string str("Chromosome ");
        str += chr;
        str += " was not found in parsed negative reads.";
        throw ChrNotFound(str);
    }
}

pritr Reads::neg_reads::end() const {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();
    return _reads._negReads.end();
}

ritr Reads::neg_reads::begin_of(string & chr) const {
    if (!_reads._noMoreNegReads)
        _reads.neg_reads._negReadsInsertionComplete();
    if (hasReadsOnChr(chr)) {
        return _reads._negReads[chr].begin();
    } else {
        string str("Chromosome ");
        str += chr;
        str += " was not found in parsed negative reads.";
        throw ChrNotFound(str);
    }
}

void Reads::neg_reads::insertRead(string & chr, uint32_t & read) {
    _reads._insertRead(chr, read, _reads._negReads);
}

uint64_t Reads::neg_reads::size() const {
    uint64_t size = 0;
    pair<string, vector<uint32_t> > reads;
    foreach(reads, _reads._negReads) {
        size += (uint64_t) (reads.second.size());
    }
    return size;
}

void Reads::neg_reads::_negReadsInsertionComplete() {
    _reads._noMoreNegReads = true;
    _reads.extractChrs(_reads._negReads, _reads._negChrs);
    LOG_DEBUG1("Neg reads insertion complete.");
    _reads.sortReadsOnEachChromosome(_reads._negReads);
    _reads._pos_is_sorted = true;

}

void Reads::neg_reads::getReads(string& chr, uint32_t start, uint32_t end,
        vector<uint32_t>::iterator& rstart, vector<uint32_t>::iterator& rend) {
    vector<uint32_t>::iterator preadsstart, preadsend;
    rt_assert_gt(end, start)
    if (hasReadsOnChr(chr)) {
        preadsstart = begin_of(chr);
        preadsend = end_of(chr);
        rstart = lower_bound(preadsstart, preadsend, start);
        rend = upper_bound(rstart, preadsend, end);
    } else {
        string str("Chromosome ");
        str += chr;
        str += " contains no negative reads.";
        throw ChrNotFound(str);
    }
}

void Reads::remove(string & chr) {
    this->pos_reads.remove(chr);
    this->neg_reads.remove(chr);
}

/*  Normal result:
 * Every chr must have reads in both the pos and neg strands,
 *  otherwise it will be removed
 *
 * Extreme result:
 * The reads are empty if reads reside in only one strand
 *
 */
void Reads::removeUnequalChrs() {
    if (!_noMoreNegReads)
        neg_reads._negReadsInsertionComplete();
    if (!_noMorePosReads)
        pos_reads._posReadsInsertionComplete();

    vector<string>::iterator it;
    vector<string> chrs_to_remove(
            neg_reads.chrs().size() + pos_reads.chrs().size());
    // todo: test this module
    if (neg_reads.chrs().size() < 1 && pos_reads.chrs().size() > 1) {
        chrs_to_remove = pos_reads.chrs();
        foreach(string chr, chrs_to_remove) {
            cerr << "Warning: " << chr
                    << " only contains reads in the positive strand."
                    << " The chromosome is removed.\n";
            pos_reads.remove(chr);
        }
        return;
    }
    if (pos_reads.chrs().size() < 1 && neg_reads.chrs().size() > 1) {
        chrs_to_remove = neg_reads.chrs();
        foreach(string chr, chrs_to_remove) {
            cerr << "Warning: " << chr
                    << " only contains reads in the negative strand."
                    << " The chromosome is removed.\n";
            neg_reads.remove(chr);
        }
        return;
    }
    //I dont trust the implementation of set_difference
    //so I have to make sure .begin() and .end() is not 0x0
    it = std::set_difference(_posChrs.begin(), _posChrs.end(), _negChrs.begin(),
            _negChrs.end(), chrs_to_remove.begin());
    chrs_to_remove.resize(it - chrs_to_remove.begin());
    foreach(string chr, chrs_to_remove) {
        cerr << "Warning: " << chr
                << " only contains reads in the positive strand."
                << " The chromosome is removed.\n";
        pos_reads.remove(chr);
    }
    vector<string> chrs_to_remove2(
            neg_reads.chrs().size() + pos_reads.chrs().size());

    it = std::set_difference(_negChrs.begin(), _negChrs.end(), _posChrs.begin(),
            _posChrs.end(), chrs_to_remove2.begin());
    chrs_to_remove2.resize(it - chrs_to_remove2.begin());
    foreach(string chr, chrs_to_remove2) {
        cerr << "Warning: " << chr
                << " only contains reads in the negative strand."
                << " The chromosome is removed.\n";
        neg_reads.remove(chr);
    }
}

