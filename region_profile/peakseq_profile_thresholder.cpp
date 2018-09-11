/*
 * peakseq_profile_thresholder.h
 *
 *  Created on: Jul 27, 2011
 *      Author: xin
 */

#include "peakseq_profile_thresholder.h"
#include "utils/logger.h"
using namespace std;
typedef pair<uint32_t, uint32_t>  profile_pos_t;

/*
 * WIll not work if the min(vec) == cutoff
 */
void peakseq_profile_thresholder::threshold(vector<profile_pos_t>& vec,
                                            uint32_t cutoff,
                                            uint32_t mergeDistance,
                                            uint32_t offset,
                                            vector<peak_pos_t>& result) {
    uint32_t peakStart = 0, peakEnd = 0, i;
    bool peakUnfinished = false;
    peak_pos_t region;

    LOG_DEBUG1("IN THRESHOLDER with parameters:"
    <<"\n size of input vec: "<<vec.size()
    <<"\n cutoff: "<<cutoff
    <<"\n offset:"<<offset
    );

    for (i = 0; i < vec.size(); i++) {
        if (vec[i].second >= cutoff && !peakUnfinished) {
            peakStart = vec[i].first;
            peakUnfinished = true;
        } else if (vec[i].second < cutoff && peakUnfinished) {
            peakUnfinished = false;
            peakEnd = vec[i].first;
            if (result.size() > 0) {
                if (peakStart - mergeDistance <= (*(result.end() - 1)).second) {
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

    LOG_DEBUG1("Total peaks found by thersholder: "<<result.size());
}
