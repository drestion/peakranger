/*
 * point_reads.cpp
 *
 *  Created on: Jul 26, 2011
 *      Author: xin
 */

#include "point_reads.h"
using namespace std;
namespace {
bool sort_comparator(pair<uint32_t, int8_t> p1, pair<uint32_t, int8_t> p2) {
    return p1.first < p2.first;
}
}
point_reads::point_reads() {

}

point_reads::~point_reads() {

}
void point_reads::set_reads(vector<uint32_t>& reads, uint32_t extension,
        vector<pos_t>& result) {
    uint32_t a;
    uint32_t b;
    uint32_t read;

    vector<uint32_t>::iterator readsStart, readsEnd;
    readsStart = reads.begin();
    readsEnd = reads.end();
    LOG_DEBUG1("start mapping reads");
    while (readsStart != readsEnd) {
        read = *readsStart++;
        a = read;
        b = read + extension;
        result.push_back(pos_t(a, 1));
        result.push_back(pos_t(b, -1));
    }

    if (result.size() == 0) {
        return;
    }
    //sort based on their locations;
    sort(result.begin(), result.end(), sort_comparator);

}

void point_reads::set_reads(uint32_t start, uint32_t end, uint32_t readlength,
        uint32_t readextlength, vector<uint32_t>::iterator readsStart,
        vector<uint32_t>::iterator readsEnd,
        vector<uint32_t>::iterator nreadsStart,
        vector<uint32_t>::iterator nreadsEnd, vector<pos_t>& result) {

    uint32_t a;
    uint32_t b;
    uint32_t arrayStart;
    uint32_t arrayEnd;
    uint32_t read;

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

                result.push_back(pair<uint32_t, int8_t>(start, 1));
                result.push_back(pair<uint32_t, int8_t>(b, -1));
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {

                result.push_back(pair<uint32_t, int8_t>(a, 1));
                result.push_back(pair<uint32_t, int8_t>(b, -1));
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {

                result.push_back(pair<uint32_t, int8_t>(a, 1));
                result.push_back(pair<uint32_t, int8_t>(end, 1));
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
        if (read + readlength < readextlength)
            continue;
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

                result.push_back(pair<uint32_t, int8_t>(start, 1));
                result.push_back(pair<uint32_t, int8_t>(b, -1));
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {

                result.push_back(pair<uint32_t, int8_t>(a, 1));
                result.push_back(pair<uint32_t, int8_t>(b, -1));
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {

                result.push_back(pair<uint32_t, int8_t>(a, 1));
                result.push_back(pair<uint32_t, int8_t>(end, 1));
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
    if (result.size() == 0) {
        return;
    }
    //sort based on their locations;
    sort(result.begin(), result.end(), sort_comparator);

}

void point_reads::set_reads(uint32_t readlength, uint32_t readextlength,
        vector<uint32_t>::iterator readsStart,
        vector<uint32_t>::iterator readsEnd, vector<pos_t>& result) {

    uint32_t a;
    uint32_t b;
    uint32_t read;

    LOG_DEBUG1("start mapping pos reads");
    while (readsStart != readsEnd) {

        read = *readsStart;
        readsStart++;
        a = read;
        b = read + readextlength;

        result.push_back(pair<uint32_t, int8_t>(a, 1));
        result.push_back(pair<uint32_t, int8_t>(b, -1));

    }

    //sort based on their locations;
    sort(result.begin(), result.end(), sort_comparator);

}

