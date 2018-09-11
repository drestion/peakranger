/*
 * PairEndedReads.h
 *
 *  Created on: May 5, 2012
 *      Author: xin
 */

#ifndef PAIRENDEDREADS_H_
#define PAIRENDEDREADS_H_
#include "ReadPair.h"
#include "common/stl_header.h"
#include "common/boost_header.h"
namespace reads {

template<typename R>
class PairEndedReads {
public:
    PairEndedReads() :
            mData() {
    }

    virtual ~PairEndedReads() {
    }

    template<typename Itr>
    void addReadPairs(Itr l, Itr r);
    void addReadPair(const ReadPair<R>& rp);

    void getReadPairs(const char* chr, std::vector<ReadPair<R> >& result);

    typename std::vector<ReadPair<R> >::iterator beginOf(
            const std::string& chr);

    typename std::vector<ReadPair<R> >::iterator endOf(const std::string& chr);

    std::vector<std::string> getChrs() const;
    size_t size();
private:

    std::map<std::string, std::vector<ReadPair<R> > > mData;

};

template<typename R>
void PairEndedReads<R>::addReadPair(const ReadPair<R>& rp) {
    mData[rp.r1().getChr()].push_back(rp);
}

template<typename R>
void PairEndedReads<R>::getReadPairs(const char* chr,
        std::vector<ReadPair<R> >& result) {
    result = mData[std::string(chr)];
}

template<typename R>
template<typename Itr>
void PairEndedReads<R>::addReadPairs(Itr l, Itr r) {
    while (l != r) {
        addReadPair(*l++);
    }
}

template<typename R>
typename std::vector<ReadPair<R> >::iterator PairEndedReads<R>::beginOf(
        const std::string& chr) {
    return mData[chr].begin();
}

template<typename R>
typename std::vector<ReadPair<R> >::iterator PairEndedReads<R>::endOf(
        const std::string& chr) {
    return mData[chr].end();
}
}

template<typename R>
inline std::vector<std::string> reads::PairEndedReads<R>::getChrs() const {
    std::vector<std::string> res;
    typedef std::map<std::string, std::vector<ReadPair<R> > > dataType;
    typename dataType::const_iterator it;

    for (it = mData.begin(); it != mData.end(); ++it) {
        res.push_back((*it).first);
    }
    return res;
}

template<typename R>
inline size_t reads::PairEndedReads<R>::size() {
    size_t res = 0;
    typedef std::map<std::string, std::vector<ReadPair<R> > > dataType;
    typename dataType::const_iterator it;

    for (it = mData.begin(); it != mData.end(); ++it) {
        res+= it->second.size();
    }
    return res;
}

/* namespace reads */
#endif /* PAIRENDEDREADS_H_ */
