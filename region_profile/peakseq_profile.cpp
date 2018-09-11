/*
 * peakseq_profile.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: xin
 */

#include "peakseq_profile.h"
#include "utils/logger.h"
using namespace std;
namespace {
bool sort_comparator(pair<uint32_t, uint32_t> p1,
                     pair<uint32_t, uint32_t> p2) {
    return p1.first < p2.first;
}

bool sort_comparator_int8(pair<uint32_t, int8_t> p1,
                          pair<uint32_t, int8_t> p2) {
    return p1.first < p2.first;
}
} /*Namespace*/

peakseq_profile::peakseq_profile() {
    // TODO Auto-generated constructor stub

}

peakseq_profile::~peakseq_profile() {
    // TODO Auto-generated destructor stub
}

void peakseq_profile::dump(vector<profile_pos_t>& profile,
                           const char* filename) {
    ofstream ofs(filename);
    for (size_t i = 0; i < profile.size(); i++) {
        ofs << profile[i].first << "\t" << profile[i].second << "\n";
    }
    ofs.close();
}

void peakseq_profile::get_profile_of_reads(vector<uint32_t> & reads,
                                           uint32_t extension,
                                           vector<profile_pos_t> & result) {
    uint32_t a;
    uint32_t b;
    uint32_t read;

    typedef vector<pair<uint32_t, int8_t> > reads_count_t;
    reads_count_t reads_count;

    vector<uint32_t>::iterator readsStart, readsEnd;
    readsStart = reads.begin();
    readsEnd = reads.end();
    LOG_DEBUG1("start mapping reads");
    while (readsStart != readsEnd) {
        read = *readsStart++;
        a = read;
        b = read + extension;
        reads_count.push_back(pair<uint32_t, int8_t> (a,
                                                      1));
        reads_count.push_back(pair<uint32_t, int8_t> (b,
                                                      -1));
    }

    if (reads_count.size() == 0) {
        return;
    }
    //sort based on their locations;
    sort(reads_count.begin(),
         reads_count.end(),
         sort_comparator_int8);

    uint32_t pos;
    uint32_t ppos = (reads_count.begin())->first;
    int32_t count = (int32_t) (reads_count.begin())->second;
    LOG_DEBUG1("start building sgr");
    reads_count_t::iterator rit = reads_count.begin();
    for (; rit != reads_count.end(); rit++) {
        pos = rit->first;
        LOG_DEBUG4("ps score:" << (uint16_t) (rit->second) << "\tps pos:"
        << pos);
        if (pos == ppos) {
            count += (int32_t) (rit->second);
            continue;
        }
        if (ppos > 0 && count > 0) {
            result.push_back(profile_pos_t(ppos,
                                           count));
        }
        count += (int32_t) (rit->second);
        ppos = pos;
    }

    if (ppos > 0 && count > 0) {
        result.push_back(profile_pos_t(ppos,
                                       count));
    }

    sort(result.begin(),
         result.end(),
         sort_comparator);

    LOG_DEBUG1("FINISHED building sgr, total profile points: "<<result.size());

}

void peakseq_profile::get_profile_of_reads(uint32_t start,
                                           uint32_t end,
                                           uint32_t readlength,
                                           uint32_t readextlength,
                                           vector<uint32_t>::iterator readsStart,
                                           vector<uint32_t>::iterator readsEnd,
                                           vector<uint32_t>::iterator nreadsStart,
                                           vector<uint32_t>::iterator nreadsEnd,
                                           vector<profile_pos_t> & result) {
    uint32_t a;
    uint32_t b;
    uint32_t read;
    typedef vector<pair<uint32_t, int8_t> > reads_count_t;
    reads_count_t reads_count;

    bool inRange = false;
    LOG_DEBUG1("start mapping pos reads");
    while (readsStart != readsEnd) {

        read = *readsStart;
        readsStart++;
        a = read;
        b = read + readextlength;

        if (a < end && b > start) {
            inRange = true;
            /*
             *     |-------|
             *   |---|
             */
            if (a <= start && b <= end) {

                reads_count.push_back(pair<uint32_t, int8_t> (start,
                                                              1));
                reads_count.push_back(pair<uint32_t, int8_t> (b,
                                                              -1));
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {

                reads_count.push_back(pair<uint32_t, int8_t> (a,
                                                              1));
                reads_count.push_back(pair<uint32_t, int8_t> (b,
                                                              -1));
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {

                reads_count.push_back(pair<uint32_t, int8_t> (a,
                                                              1));
                reads_count.push_back(pair<uint32_t, int8_t> (end,
                                                              1));
            }
            /*
             *     |-------|
             *   |-----------|
             */
            if (a < start && b > end) {

            }

        } else if (inRange) {
            break;
        }
    }

    LOG_DEBUG1("start mapping neg reads");

    while (nreadsStart != nreadsEnd) {

        read = *nreadsStart;
        nreadsStart++;
        if (read + readlength < readextlength) continue;
        a = read + readlength - readextlength;
        b = read;

        if (a < end && b > start) {
            inRange = true;
            /*
             *     |-------|
             *   |---|
             */
            if (a <= start && b <= end) {

                reads_count.push_back(pair<uint32_t, int8_t> (start,
                                                              1));
                reads_count.push_back(pair<uint32_t, int8_t> (b,
                                                              -1));
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {

                reads_count.push_back(pair<uint32_t, int8_t> (a,
                                                              1));
                reads_count.push_back(pair<uint32_t, int8_t> (b,
                                                              -1));
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {

                reads_count.push_back(pair<uint32_t, int8_t> (a,
                                                              1));
                reads_count.push_back(pair<uint32_t, int8_t> (end,
                                                              1));
            }
            /*
             *     |-------|
             *   |-----------|
             */
            if (a < start && b > end) {

            }

        } else if (inRange) {
            break;
        }
    }
    if (reads_count.size() == 0) {
        return;
    }
    //sort based on their locations;
    sort(reads_count.begin(),
         reads_count.end(),
         sort_comparator_int8);

    uint32_t pos;
    uint32_t ppos = (reads_count.begin())->first;
    int32_t count = (int32_t) (reads_count.begin())->second;
    LOG_DEBUG1("start building sgr");
    reads_count_t::iterator rit = reads_count.begin();
    for (; rit != reads_count.end(); rit++) {
        pos = rit->first;
        LOG_DEBUG4("ps score:" << (uint16_t) (rit->second) << "\tps pos:"
        << pos);
        if (pos == ppos) {
            count += (int32_t) (rit->second);
            continue;
        }
        if (ppos > 0 && count > 0) {
            result.push_back(profile_pos_t(ppos,
                                           count));
        }
        count += (int32_t) (rit->second);
        ppos = pos;
    }

    if (ppos > 0 && count > 0) {
        result.push_back(profile_pos_t(ppos,
                                       count));
    }

    sort(result.begin(),
         result.end(),
         sort_comparator);

    LOG_DEBUG1("FINISHED building sgr, total profile points: "<<result.size());

}

