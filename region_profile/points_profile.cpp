/*
 * points_profile.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: xin
 */

#include "points_profile.h"
#include "utils/logger.h"
#include "utils/assert_helpers.h"
using namespace std;

namespace {
bool sort_comparator(pair<uint32_t, int8_t> p1,
                     pair<uint32_t, int8_t> p2) {
    return p1.first < p2.first;
}
}
points_profile::points_profile() {
    // TODO Auto-generated constructor stub

}

points_profile::~points_profile() {
    // TODO Auto-generated destructor stub
}
/*
 * Reads must be sorted before calling this method
 */
void points_profile::get_profile_of_reads(vector<uint32_t>& reads,
                                          uint32_t extension,
                                          vector<uint16_t>& result) {
    uint32_t a;
    uint32_t b;
    uint32_t arrayStart = *(reads.begin());
    uint32_t arrayEnd = *(reads.end() - 1) + extension+1;
    uint32_t arrayLength = arrayEnd - arrayStart + 1;
    uint32_t read;

    typedef vector<pair<uint32_t, int8_t> > reads_count_t;
    reads_count_t reads_count;
    result.resize(arrayLength);
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
         sort_comparator);

    uint32_t pos;
    uint32_t ppos = (reads_count.begin())->first;
    uint16_t count = (reads_count.begin())->second;
    LOG_DEBUG1("start building sgr");
    reads_count_t::iterator rit = reads_count.begin();
    for (; rit != reads_count.end(); rit++) {
        pos = rit->first;
        LOG_DEBUG4("ps score:" << (uint16_t) (rit->second)<<"\tps pos:" << pos);
        if (pos == ppos) {
            count += (uint16_t) (rit->second);
            continue;
        }
        if (ppos > 0 && count > 0) {
            result.at(ppos - arrayStart) = (uint16_t) count;
        }
        count += (uint16_t) (rit->second);
        ppos = pos;
    }

    LOG_DEBUG1("FINISHED building sgr");
    if (ppos > 0 && count > 0) {
        result.at(ppos - arrayStart) = (uint16_t) count;
    }
}
/*
 * Reads must be sorted before calling this method
 */
void points_profile::get_profile_of_reads(uint32_t start,
                                          uint32_t end,
                                          uint32_t readlength,
                                          uint32_t readextlength,
                                          vector<uint32_t>::iterator readsStart,
                                          vector<uint32_t>::iterator readsEnd,
                                          vector<uint32_t>::iterator nreadsStart,
                                          vector<uint32_t>::iterator nreadsEnd,
                                          vector<uint16_t>& result) {

    uint32_t a;
    uint32_t b;
    uint32_t arrayStart;
    uint32_t arrayEnd;
    uint32_t arrayLength = end - start + 1;
    uint32_t read;
    typedef vector<pair<uint32_t, int8_t> > reads_count_t;
    reads_count_t reads_count;
    result.resize(arrayLength);

    bool inRange = false;
    LOG_DEBUG1("start mapping pos reads");
    while (readsStart != readsEnd) {

        read = *readsStart;
        readsStart++;
        a = read;
        b = read + readextlength;
        arrayStart = 0;
        arrayEnd = 0;

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
        arrayStart = 0;
        arrayEnd = 0;

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
         sort_comparator);

    uint32_t pos;
    uint32_t ppos = (reads_count.begin())->first;
    uint16_t count = (reads_count.begin())->second;
    LOG_DEBUG1("start building sgr");
    reads_count_t::iterator rit = reads_count.begin();
    for (; rit != reads_count.end(); rit++) {
        pos = rit->first;
        LOG_DEBUG4("ps score:" << (uint16_t) (rit->second)<<"\tps pos:" << pos);
        if (pos == ppos) {
            count += (uint16_t) (rit->second);
            continue;
        }
        if (ppos > 0 && count > 0) {
            //            os << ppos << "\t" << count << "\n";
            result.at(ppos - start) = (uint16_t) count;
        }
        count += (uint16_t) (rit->second);
        ppos = pos;
    }

    LOG_DEBUG1("FINISHED building sgr");
    if (ppos > 0 && count > 0) {
        //        os << ppos << "\t" << count << "\n";
        result.at(ppos - start) = (uint16_t) count;
    }
}
