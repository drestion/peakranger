/*
 * Strand.h
 *
 *  Created on: May 3, 2012
 *      Author: xin
 */

#ifndef STRAND_H_
#define STRAND_H_

#include <ostream>

namespace reads {

class Strand {
    friend std::ostream& operator<<(std::ostream& os, const Strand& rhs) {
        if (rhs.isPos()) {
            os << "+";
        } else {
            os << "-";
        }
        return os;
    }
public:
    Strand();
    Strand(bool pos);
    Strand(const Strand& str);
    virtual ~Strand();
    void operator=(bool dir);
    bool operator==(const Strand& rhs) const;
    bool isPos() const;
    bool isNeg() const;
    void toPos();
    void toNeg();

private:
    bool mDir;
};

} /* namespace reads */
#endif /* STRAND_H_ */
