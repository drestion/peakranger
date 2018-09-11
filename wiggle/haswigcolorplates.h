/*
 * haswigcolorplates.h
 *
 *  Created on: Jan 26, 2012
 *      Author: xfeng
 */

#ifndef HASWIGCOLORPLATES_H_
#define HASWIGCOLORPLATES_H_
#include <ostream>
#include <vector>
#include <string>
#include <stdint.h>

class has_wig_color_plates{
public:
    has_wig_color_plates();
    virtual ~has_wig_color_plates();

protected:
    std::vector<uint32_t> _red;
    std::vector<uint32_t> _brown;
    std::vector<uint32_t> _green;

};

#endif /* HASWIGCOLORPLATES_H_ */
