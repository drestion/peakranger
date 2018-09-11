/*
 * reads.h
 *
 *  Created on: Apr 28, 2011
 *      Author: xin
 */

#ifndef READS_H_
#define READS_H_

#include "utils/exceptions.h"
#include "utils/logger.h"
#include <vector>
#include <utility>
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


typedef std::vector<uint32_t> reads_vec;
typedef std::map<std::string, reads_vec> reads_t;
typedef std::map<std::string, reads_vec>::iterator pritr;
typedef std::vector<uint32_t>::iterator ritr;

/*
 * Reads are categorized based on strands and chromosomes.
 * No assumptions are posed toward strands and chromosomes.
 * It is possible that a chr contains only pos or neg reads.
 *
 *
 *
 */
class Reads {
public:
    Reads() :
            pos_reads(*this), neg_reads(*this), _noMorePosReads(false), _noMoreNegReads(
                    false), _pos_is_sorted(false), _neg_is_sorted(false) {
    }

    struct pos_reads {
    public:
        friend class Reads;

        /*
         * access to all reads;
         */
        pritr begin() const;
        pritr end() const;

        /*
         * access to reads of a chr
         */
        ritr begin_of(std::string& chr) const;
        ritr end_of(std::string& chr) const;

        /*
         * Access to reads in a region
         */
        void getReads(std::string& chr, uint32_t start, uint32_t end,
                std::vector<uint32_t>::iterator& rstart,
                std::vector<uint32_t>::iterator& rend);
        /*
         * reads insertion
         */
        void insertRead(std::string& chr, uint32_t& read);
        /*
         * properties
         */
        uint64_t size() const;
        std::vector<std::string> chrs() const;
        bool hasReadsOnChr(std::string& chr) const;

    private:
        inline void _posReadsInsertionComplete();
        Reads& _reads;
        pos_reads(Reads& reads) :
                _reads(reads) {
        }
        /*
         * modification
         */
        void remove(std::string& chr);

    } pos_reads;

    struct neg_reads {
    public:
        friend class Reads;
        typedef std::map<std::string, std::vector<uint32_t> >::iterator pritr;
        /*
         * access to all reads;
         */
        pritr begin() const;
        pritr end() const;

        /*
         * access to reads of a chr
         */
        ritr begin_of(std::string& chr) const;
        ritr end_of(std::string& chr) const;

        /*
         * Access to reads in a region
         */
        void getReads(std::string& chr, uint32_t start, uint32_t end,
                std::vector<uint32_t>::iterator& rstart,
                std::vector<uint32_t>::iterator& rend);

        /*
         * reads insertion
         */
        void insertRead(std::string& chr, uint32_t& read);
        /*
         * properties
         */
        uint64_t size() const;
        std::vector<std::string> chrs() const;
        bool hasReadsOnChr(std::string& chr) const;

    private:
        inline void _negReadsInsertionComplete();
        neg_reads(Reads& reads) :
                _reads(reads) {
        }
        /*
         * modification
         */
        void remove(std::string& chr);
        Reads& _reads;

    } neg_reads;

    /*
     * return the number of pos and neg reads
     */
    uint64_t size() const;

    /*
     * Properties
     */
    uint32_t getReadlength() const;
    void setReadlength(uint32_t _readlength);

    /*
     * modification
     */
    void remove(std::string& chr);

    /*
     * Quality control
     */
    void removeUnequalChrs();

protected:
    void _insertRead(std::string& chr, uint32_t read, reads_t& reads);

    void sortReadsOnEachChromosome(reads_t& reads) {
        reads_t::iterator itr = reads.begin();
        for (; itr != reads.end(); itr++) {
            LOG_DEBUG3("Sorted "<<itr->first);
            sort(itr->second.begin(), itr->second.end());
        }
    }

    void extractChrs(reads_t& reads, std::vector<std::string>& chrs) {
        std::pair<std::string, reads_vec> r;
        foreach(r, reads) {
            chrs.push_back(r.first);
        }
        sort(chrs.begin(), chrs.end());
    }

protected:
    std::map<std::string, std::vector<uint32_t> > _posReads;
    std::map<std::string, std::vector<uint32_t> > _negReads;
    std::vector<std::string> _negChrs;
    std::vector<std::string> _posChrs;
    uint32_t _readlength;
    bool _noMorePosReads;
    bool _noMoreNegReads;
    bool _pos_is_sorted;
    bool _neg_is_sorted;

};
#endif /* READS_H_ */
