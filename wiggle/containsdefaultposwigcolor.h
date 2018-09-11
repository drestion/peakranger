/*
 * containsdefaultposwigcolor.h
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#ifndef CONTAINSDEFAULTPOSWIGCOLOR_H_
#define CONTAINSDEFAULTPOSWIGCOLOR_H_
#include <ostream>
#include <vector>
#include <string>
#include <stdint.h>


class contains_default_pos_wig_color{
public:
    contains_default_pos_wig_color();
    virtual ~contains_default_pos_wig_color();
    std::vector<uint32_t> getPosRgb() const;
    void setPosRgb(std::vector<uint32_t> _posRGB);

private:
    std::vector<uint32_t> _posRGB;
};

#endif /* CONTAINSDEFAULTPOSWIGCOLOR_H_ */
