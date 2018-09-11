/*
 * readstools.cpp
 *
 *  Created on: Sep 22, 2011
 *      Author: xfeng
 */

#include "readstools.h"

#define MAX_CHR_NUM 100000
using namespace std;
reads_tools::reads_tools()
{
    
}

reads_tools::~reads_tools()
{

}

void reads_tools::verify_and_correct_Reads_both_strands(Reads & treads,
                                                        Reads & creads) {
    vector<string>::iterator it;

    treads.removeUnequalChrs();
    creads.removeUnequalChrs();
    vector < string > treads_chrs = treads.pos_reads.chrs(); //sorted
    vector < string > creads_chrs = creads.pos_reads.chrs(); //sorted
    //Must initialize it before using it.
    vector < string > chrs_to_remove(treads_chrs.size() + creads_chrs.size());
    if (chrs_to_remove.size() < 1) {
        return;
    }

    it = std::set_difference(treads_chrs.begin(),
                             treads_chrs.end(),
                             creads_chrs.begin(),
                             creads_chrs.end(),
                             chrs_to_remove.begin());
    chrs_to_remove.resize(it - chrs_to_remove.begin());
    /*
     * Use set_symmertric_difference is you want them at one time
     */foreach(string chr, chrs_to_remove)
    {

        cerr << "Warning: "
        << "No reads were found in " << chr
        << " of the control dataset. The chromosome is removed.\n";

        treads.remove(chr);
    }
    vector<string>::iterator it2;
    //Must initialize it before using it.
    vector < string > chrs_to_remove2(treads_chrs.size() + creads_chrs.size());
    it2 = std::set_difference(creads_chrs.begin(),
                              creads_chrs.end(),
                              treads_chrs.begin(),
                              treads_chrs.end(),
                              chrs_to_remove2.begin());
    chrs_to_remove2.resize(it2 - chrs_to_remove2.begin());
    foreach(string chr, chrs_to_remove2)
    {
        cerr << "Warning: "
        << "No reads were found in " << chr
        << " of the treatment dataset. The chromosome is removed.\n";
        creads.remove(chr);
    }
}

void reads_tools::get_merged_chrs_for_both_strands(Reads & reads,
                                                   vector<string> & mergedchrs)
                                                   {
    if (reads.size() < 1) {
        return;
    }
    vector < string > pchrs = reads.pos_reads.chrs();
    vector < string > nchrs = reads.neg_reads.chrs();
    if (mergedchrs.size() < pchrs.size() + nchrs.size()) {
        mergedchrs.resize(pchrs.size() + nchrs.size());
    }

    vector<string>::iterator chrit;
    if (pchrs.size() > 1)
    {
        sort(pchrs.begin(),
             pchrs.end());
    }
    if (nchrs.size() > 1) {
        sort(nchrs.begin(),
             nchrs.end());
    }
    chrit = merge(pchrs.begin(),
                  pchrs.end(),
                  nchrs.begin(),
                  nchrs.end(),
                  mergedchrs.begin());
    mergedchrs.resize(chrit - mergedchrs.begin());
    chrit = unique(mergedchrs.begin(),
                   mergedchrs.end());
    mergedchrs.resize(chrit - mergedchrs.begin());
}

void reads_tools::_insert_random_reads(Reads & reads,
                                       string chr,
                                       Reads & result)
                                       {
    vector<uint32_t> _rds;
    uint32_t _rs = reads.pos_reads.size();

    uint32_t _l = (reads.pos_reads.begin_of(chr))[0];
    uint32_t _h = (reads.pos_reads.end_of(chr) - 1)[0];
    _rs = reads.pos_reads.end_of(chr) - reads.pos_reads.begin_of(chr);
    if (_rs < 1) {
        return;
    }
    distributions::random_int_std(_l,
                                  _h,
                                  _rs,
                                  _rds);
    foreach(uint32_t _rd, _rds)
    {
        result.pos_reads.insertRead(chr,
                                    _rd);
    }
}

void reads_tools::generate_random_reads_based_on_reads(Reads & reads,
                                                       Reads & result)
                                                       {
    vector < string > chrs;
    vector<uint32_t> _rds;
    //todo: shall this be an option?
    reads.removeUnequalChrs();

    get_merged_chrs_for_both_strands(reads,
                                     chrs);

    foreach(string chr, chrs)
    {

        _insert_random_reads(reads,
                             chr,
                             result);

        _rds.clear();
        uint32_t _rs = reads.neg_reads.size();
        uint32_t _l = (reads.neg_reads.begin_of(chr))[0];
        uint32_t _h = (reads.neg_reads.end_of(chr) - 1)[0];
        _rs = reads.neg_reads.end_of(chr) - reads.neg_reads.begin_of(chr);
        if (_rs < 1) continue;
        distributions::random_int_std(_l,
                                      _h,
                                      _rs,
                                      _rds);
        foreach(uint32_t _rd, _rds)
        {
            result.neg_reads.insertRead(chr,
                                        _rd);
        }
    }

    result.setReadlength(reads.getReadlength());
    result.pos_reads.begin();
    result.neg_reads.begin();
}

size_t reads_tools::chromSize(Reads & treads,
                               string& chr)
                              {
    return *(treads.pos_reads.end_of(chr) - 1)
    > *(treads.neg_reads.end_of(chr) - 1) ?
    *(treads.pos_reads.end_of(chr) - 1)
    + treads.getReadlength() + 1000 :
    *(treads.neg_reads.end_of(chr) - 1)
    + treads.getReadlength() + 1000;
}

