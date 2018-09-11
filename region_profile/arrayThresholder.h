/*
 * arrayThresholder.h
 *
 *  Created on: May 3, 2011
 *      Author: xin
 */

#ifndef ARRAYTHRESHOLDER_H_
#define ARRAYTHRESHOLDER_H_

#include <vector>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include <stdint.h>
#include "utils/exceptions.h"

template<class T>
class arrayThresholder {
public:

    arrayThresholder() {
    }
    ~arrayThresholder() {
    }

public:
    /*
     * result will be reset
     */
    static void threshold(std::vector<T>& vec, T cutoff, uint32_t mergeDistance,
            std::vector<
                    std::pair<typename std::vector<T>::iterator,
                            typename std::vector<T>::iterator> >& result) {
        /*
         * result.size() must be 0
         */
        result.clear();

        bool peakUnfinished = false;
        typename std::vector<T>::iterator peakStart, peakEnd, i;
        std::pair<typename std::vector<T>::iterator,
                typename std::vector<T>::iterator> region;

        for (i = vec.begin(); i != vec.end(); i++) {
            if (*i > cutoff && !peakUnfinished) {
                peakStart = i;
                peakUnfinished = true;
            } else if (*i < cutoff && peakUnfinished) {
                peakUnfinished = false;
                peakEnd = i;
                if (result.size() > 0) {
                    if (peakStart - mergeDistance
                            <= (*(result.end() - 1)).second) {
                        (*(result.end() - 1)).second = peakEnd;
                        continue;
                    }
                }
                region.first = peakStart;
                region.second = peakEnd;
                result.push_back(region);
            }
        }

        if (peakUnfinished) {
            region.first = peakStart;
            region.second = vec.end() - 1;
            result.push_back(region);
        }

    }

    /*
     * result will be reset
     */
    static void threshold(std::vector<T>& vec, T cutoff, uint32_t mergeDistance,
            uint32_t offset, std::vector<std::pair<uint32_t, uint32_t> >& result) {
        /*
         * result.size() must be 0
         */
        result.clear();

        bool peakUnfinished = false;
        uint32_t peakStart = 0;
        uint32_t peakEnd = 0;
        uint32_t i;
        std::pair<uint32_t, uint32_t> region;

        for (i = 0; i < vec.size(); i++) {
            if (vec[i] > cutoff && !peakUnfinished) {
                peakStart = i;
                peakUnfinished = true;
            } else if (vec[i] < cutoff && peakUnfinished) {
                peakUnfinished = false;
                peakEnd = i;
                if (result.size() > 0) {
                    if (peakStart - mergeDistance
                            <= (*(result.end() - 1)).second) {
                        (*(result.end() - 1)).second = peakEnd;
                        continue;
                    }
                }
                region.first = peakStart;
                region.second = peakEnd;
                result.push_back(region);
            }
        }

        if (peakUnfinished) {
            region.first = peakStart;
            region.second = vec.size() - 1;
            result.push_back(region);
        }

        for (i = 0; i < result.size(); i++) {
            result[i].first += offset;
            result[i].second += offset;
        }
    }

    static void threshold_mr(std::vector<T>& vec, T cutoff,
            uint32_t mergeDistance,
            std::vector<
                    std::pair<typename std::vector<T>::iterator,
                            typename std::vector<T>::iterator> >& result) {
        bool peakUnfinished = false;
        typename std::vector<T>::iterator peakStart, peakEnd, i;
        std::pair<typename std::vector<T>::iterator,
                typename std::vector<T>::iterator> region;
        for (i = vec.begin(); i != vec.end(); i++) {
            if (*i > cutoff && !peakUnfinished) {
                peakStart = i;
                peakUnfinished = true;
            } else if (*i < cutoff && peakUnfinished) {
                peakUnfinished = false;
                peakEnd = i;
                if (result.size() > 0) {
                    if (peakStart - mergeDistance
                            <= (*(result.end() - 1)).second) {
                        (*(result.end() - 1)).second = peakEnd;
                        continue;
                    }
                }
                region.first = peakStart;
                region.second = peakEnd;
                result.push_back(region);
            }
        }

        if (peakUnfinished) {
            region.first = peakStart;
            region.second = vec.end() - 1;
            result.push_back(region);
        }

    }

};

#endif /* ARRAYTHRESHOLDER_H_ */
