/*
 * readstools.h
 *
 *  Created on: Sep 22, 2011
 *      Author: xfeng
 */

#ifndef READSTOOLS_H_
#define READSTOOLS_H_
#include "reads.h"
#include "utils/assert_helpers.h"
#include "utils/exceptions.h"
#include "math/distributions.h"
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include "common/stl_header.h"

#define foreach BOOST_FOREACH
/*
 * Manipulate reads
 */
class reads_tools {
public:
    reads_tools();
    virtual ~reads_tools();
    /*
     * Remove chrs that only have reads in one chr.
     * Also remove chrs that only exist in one Reads object
     */
    static void verify_and_correct_Reads_both_strands(Reads & treads,
            Reads & creads);

    /*
     * Get all chrs ever appeared in both strands
     */
    static void get_merged_chrs_for_both_strands(Reads& reads,
            std::vector<std::string>& mergedchrs);

    static void generate_random_reads_based_on_reads(Reads& reads,
            Reads& result);

    static size_t chromSize(Reads& reads, std::string& chr);

private:
    static void _insert_random_reads(Reads & reads, std::string chr, Reads & result);
};

#endif /* READSTOOLS_H_ */
