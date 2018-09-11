/*
 * TabGene.h
 *
 *  Created on: Mar 8, 2012
 *      Author: xinfeng
 */

#ifndef TABGENE_H_
#define TABGENE_H_

#include <vector>
#include <utility>
#include <string>
#include <iostream>
#include "Region.h"

namespace tab_file {

class TabGene{
    friend std::ostream& operator<<(std::ostream& os,
                                    const TabGene& tg);
    public:
    TabGene();
    TabGene(const Region& utr5,
            const Region& utr3,
            const std::vector<Region>& exons,
            const std::string& name,
            bool dir);
    virtual ~TabGene();
    bool operator<(const TabGene& rhs) const;

    Region utr5;
    Region utr3;
    std::vector<Region> exons;
    bool dir;
    std::string name;


};

bool sameID(const TabGene& lhs,
            const TabGene& rhs);

bool name_lex_lt(const TabGene& lhs,
                 const TabGene& rhs);

bool overlaps(const TabGene& lhs,
              const TabGene& rhs);
}

/* namespace tab_file */
#endif /* TABGENE_H_ */
