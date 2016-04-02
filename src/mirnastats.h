#ifndef _MIRNA_STATS_H_
#define _MIRNA_STATS_H_

#include <vector>
#include <string>

struct miRNAstats
{
    double n_stems, n_asyms, n_syms, n_bases, n_base_pairs, n_loop_bases;
    std::vector < double >stem_len, asym_len, asym_nbases, sym_len;
};                              //This structure includes information about the structure of a given miRNA sequence


miRNAstats getstats (std::string inp, bool considertruncated=false);     //For Matrix

#endif
