/*
 * wig.h
 *
 *  Created on: Jan 12, 2012
 *      Author: xfeng
 */

#ifndef WIG_H_
#define WIG_H_
#include <stdint.h>
#include <vector>
class wig{
public:
    wig();
    wig(uint32_t p,
        double s);
    uint32_t getP() const;
    double getS() const;
    void setP(uint32_t p);
    void setS(double s);
    bool operator<(const wig& w);
    static bool compa(const wig& w1,
                      const wig& w2);

protected:
    uint32_t _p;
    double _s;
};

typedef std::vector<wig> wigs;

#endif /* WIG_H_ */
