/*
 * TabGene.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: xinfeng
 */

#include "TabGene.h"

namespace tab_file {

TabGene::TabGene()
: utr5(Region(0,
              0)),
  utr3(Region(0,
              0)),
  exons(), dir(true), name("") {

}

TabGene::TabGene(const Region& utr5,
                 const Region& utr3,
                 const std::vector<Region>& exons,
                 const std::string& name,
                 bool dir)
: utr5(utr5), utr3(utr3), exons(exons), dir(dir), name(name) {
}

TabGene::~TabGene() {

}

std::ostream& operator <<(std::ostream & os,
                          const TabGene& tg)
                          {

    os << "name2:" << tg.name << ", "
    << "tx:(" << tg.utr5.getL() << ", " << tg.utr3.getR() << "), "
    << "cds:(" << tg.utr5.getR() << ", " << tg.utr3.getL() << "), "
    << "exons:(";
    for (size_t i = 0; i < tg.exons.size(); i++) {
        os << tg.exons[i].getL() << "-" << tg.exons[i].getR() << ", ";
    }
    os << "), "
    << "utr5:(" << tg.utr5 << "), "
    << "utr3:(" << tg.utr3 << ")"
    << "\n";
    return os;
}

bool TabGene::operator <(const TabGene & rhs) const
                         {
    Region tx(rhs.utr5.getL(),
              rhs.utr3.getR());
    Region txm(utr5.getL(),
               utr3.getR());
    return txm < tx;
}

bool name_lex_lt(const TabGene & lhs,
                 const TabGene & rhs)
                 {
    int ret = lhs.name.compare(rhs.name);
    return ret < 0;
}

}/* namespace tab_file */

bool tab_file::sameID(const TabGene & lhs,
                      const TabGene & rhs)
                      {
    return lhs.name == rhs.name;
}

bool tab_file::overlaps(const TabGene & lhs,
                        const TabGene & rhs)
                        {
    Region a(rhs.utr5.getL(),
             rhs.utr3.getR());
    Region b(lhs.utr5.getL(),
             lhs.utr3.getR());

    return b.overlaps(a);

}
