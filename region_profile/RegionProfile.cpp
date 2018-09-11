/*
 * RegionProfile.cpp
 *
 *  Created on: Apr 27, 2011
 *      Author: xin
 */

#include "RegionProfile.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include <iostream>
using namespace std;

/*
 * Is more like one designed exclusively for PeakSeq
 */
void region_profile::get_region_profile(uint32_t start,
                                        uint32_t end,
                                        uint32_t extension,
                                        vector<uint32_t>::iterator readsStart,
                                        vector<uint32_t>::iterator readsEnd,
                                        vector<uint16_t> & result) {
    assert_gt(end,start)
    uint32_t a;
    uint32_t b;
    uint32_t arrayStart;
    uint32_t arrayEnd;
    uint32_t arrayLength = end - start + 1;
    uint32_t read;

    result.resize(arrayLength);

    bool inRange = false;

    while (readsStart != readsEnd) {

        read = *readsStart;
        readsStart++;
        a = read;
        b = read + extension;
        arrayStart = 0;
        arrayEnd = 0;

        if (a < end && b > start) {
            inRange = true;
            /*
             *     |-------|
             *   |---|
             */
            if (a <= start && b <= end) {
                arrayStart = 0;
                arrayEnd = (b - start + 1);
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {
                arrayStart = (a - start);
                arrayEnd = (b - start + 1);
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {
                arrayStart = (a - start);
                arrayEnd = arrayLength;
            }
            /*
             *     |-------|
             *   |-----------|
             */
            if (a < start && b > end) {
                arrayStart = 0;
                arrayEnd = arrayLength;
            }

            while (arrayStart < arrayEnd) {
                ++result[arrayStart];
                arrayStart++;
            }
        } else if (inRange) {
            break;
        }
    }
}

/*
 * Is more like one designed exclusively for PeakSeq
 *
 * defers only in the type of result
 */
void region_profile::get_region_profile(uint32_t start,
                                        uint32_t end,
                                        uint32_t extension,
                                        vector<uint32_t>::iterator readsStart,
                                        vector<uint32_t>::iterator readsEnd,
                                        vector<uint64_t> & result) {

    uint32_t a;
    uint32_t b;
    uint32_t arrayStart;
    uint32_t arrayEnd;
    uint32_t arrayLength = end - start + 1;
    uint32_t read;

    result.resize(arrayLength);

    bool inRange = false;

    while (readsStart != readsEnd) {

        read = *readsStart;
        readsStart++;
        a = read;
        b = read + extension;
        arrayStart = 0;
        arrayEnd = 0;

        if (a < end && b > start) {
            inRange = true;
            /*
             *     |-------|
             *   |---|
             */
            if (a <= start && b <= end) {
                arrayStart = 0;
                arrayEnd = (b - start + 1);
            }
            /*
             *     |-------|
             *       |---|
             */
            if (a >= start && b <= end) {
                arrayStart = (a - start);
                arrayEnd = (b - start + 1);
            }
            /*
             *     |-------|
             *            |---|
             */
            if (a < end && b > end) {
                arrayStart = (a - start);
                arrayEnd = arrayLength;
            }
            /*
             *     |-------|
             *   |-----------|
             */
            if (a < start && b > end) {
                arrayStart = 0;
                arrayEnd = arrayLength;
            }

            while (arrayStart < arrayEnd) {
                ++result[arrayStart];
                arrayStart++;
            }
        } else if (inRange) {
            break;
        }
    }
}

/*
 * result will be cleared before profiling.
 *
 * reads must be sorted
 */
void region_profile::get_profile_of_reads(vector<uint32_t> & reads,
                                          uint32_t extension,
                                          vector<uint16_t> & result) {

    if (reads.begin() == reads.end()) return;
    if (reads.size() < 1) return;
    uint32_t a;
    uint32_t b;
    uint32_t arrayStart = *(reads.begin());
    uint32_t arrayEnd = *(reads.end() - 1) + extension;
    uint32_t arrayLength = arrayEnd - arrayStart + 1;
    uint32_t read;
    vector<uint32_t>::iterator it = reads.begin();
    result.resize(arrayLength);

    while (it != reads.end()) {
        read = *it++;
        a = read;
        b = read + extension;
        while (a < b) {
           if( a - *(reads.begin()) > 0)
            ++result[a - *(reads.begin())];
            a++;
        }
    }
}

