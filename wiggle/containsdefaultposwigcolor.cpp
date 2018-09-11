/*
 * containsdefaultposwigcolor.cpp
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#include "containsdefaultposwigcolor.h"
#include <ostream>
#include <vector>
#include <string>
#include <stdint.h>

using namespace std;

contains_default_pos_wig_color::contains_default_pos_wig_color()
{
    uint32_t ga[3] = { 225, 127, 0 };
    _posRGB = vector < uint32_t > (ga,
    ga + 3);
    
}

contains_default_pos_wig_color::~contains_default_pos_wig_color()
{

}

vector<uint32_t> contains_default_pos_wig_color::getPosRgb() const
{
    return _posRGB;
}

void contains_default_pos_wig_color::setPosRgb(vector<uint32_t> _posRGB)
                                               {
    this->_posRGB = _posRGB;
}

