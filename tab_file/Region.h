/*
 * Region.h
 *
 *  Created on: Mar 8, 2012
 *      Author: xinfeng
 */

#ifndef REGION_H_
#define REGION_H_

#include <cstring>
#include <iostream>
namespace tab_file {
//@to be depreciated
class Region{
    friend std::ostream& operator<<(std::ostream& os,
                                    const Region& tg);

    public:
    typedef size_t loc_t;
    Region();
    Region(loc_t l,
           loc_t r);
    Region(const Region& r);

    virtual ~Region();

    bool operator<(const Region& rhs) const;
    bool operator>(const Region& rhs) const;
    bool operator==(const Region& rhs) const;
    bool operator[](const Region& rhs) const;
    bool overlaps(const Region& rhs) const;
    Region& operator=(const Region& rhs);
    loc_t getL() const;
    void setL(loc_t l);
    void setR(loc_t r);
    loc_t getR() const;
    protected:

private:

    loc_t l;
    loc_t r;
    void swap(Region& r);
};

} /* namespace tab_file */
#endif /* REGION_H_ */
