/*
 * RegionProfile.h
 *
 *  Created on: Apr 27, 2011
 *      Author: xin
 */

#ifndef REGIONPROFILE_H_
#define REGIONPROFILE_H_

#include <stdint.h>
#include <algorithm>
#include <vector>
#include "short_reads/reads.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

class region_profile {

public:

    /*
     * It will transform the (neg) read using
     * trans(read). Then it just performs the same job
     * as other profilers. We included redudant codes
     * here but we saved memory.
     *
     * The transform function should be in the form
     *
     * trans(read_t read, uint32_t, read_length,uint32_t extension);
     *
     * The content of result will be kept for 0~(end-start+1)
     *
     * Its's the user's responsibility to properly initialize
     * the result vector
     */
    template<class neg_read_transformer>
    static void get_region_profile(uint32_t start, uint32_t end,
            uint32_t extension, uint32_t read_length,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint16_t>& result, neg_read_transformer trans) {
        assert_gt(end, start)
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
            /* transform the read */

            a = trans(read, read_length, extension);

            /* end transform */
            b = a + extension; // use the transformed read
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
     * style. used to process the pos reads
     */
    static void
    get_region_profile(uint32_t start, uint32_t end, uint32_t extension,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint16_t>& result);

    /*
     * Is more like one designed exclusively for PeakSeq
     * style. used to process the pos reads
     */
    static void
    get_region_profile(uint32_t start, uint32_t end, uint32_t extension,
            std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint64_t>& result);

    /*
     * Process both pos and neg strands.
     */
    template<class neg_read_transformer>
    static void get_region_profile(Reads& reads, std::string& chr,
            uint32_t start, uint32_t end, uint32_t extension,
            std::vector<uint16_t>& result, neg_read_transformer trans) {

        typedef std::vector<uint32_t>::iterator readitr_t;
        readitr_t readsStart, readsEnd;
        readsStart = reads.pos_reads.begin_of(chr);
        readsEnd = reads.pos_reads.end_of(chr);
        uint32_t readlength = reads.getReadlength();

        /*
         * pos strand
         */
        get_region_profile(start, end, extension, readsStart, readsEnd, result);

        readsStart = reads.neg_reads.begin_of(chr);
        readsEnd = reads.neg_reads.end_of(chr);

        /*
         * neg strand
         */
        get_region_profile(start, end, extension, readlength, readsStart,
                readsEnd, result, trans);
    }

    /*
     * Get the profile using reads on
     * both strands in the PeakSeq fasion,
     * which convert neg reads to pos reads
     * by:
     * t_read = neg_read + readlength - extension
     */
    static void get_region_profile_peakseq(Reads& reads, std::string& chr,
            uint32_t start, uint32_t end, uint32_t extension,
            std::vector<uint16_t>& result) {
        get_region_profile(reads, chr, start, end, extension, result,
                neg_read_trans);
    }

    static void get_profile_of_reads(std::vector<uint32_t>& reads,
            uint32_t extension, std::vector<uint16_t>& result);
private:
    /*
     * the basic neg read transfomer using
     * the algorithm from PeakSeq
     */
    static uint32_t neg_read_trans(uint32_t read, uint32_t readlength,
            uint32_t extension) {
        uint32_t sum = read + readlength;
        if (sum > extension) {
            return sum - extension;
        } else
            return 0;
    }

};

#endif /* REGIONPROFILE_H_ */
