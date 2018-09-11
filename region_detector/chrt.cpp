/*
 * chrt.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#include "chrt.h"

namespace ccat_aux {

chr_t::chr_t()
: chromName(""),
  chromIndex(0),
  chromSize(0),
  l1PosTags(),
  l1NegTags(),
  l2PosTags(),
  l2NegTags(),
  l1Peaks(),
  l2Peaks(0)
{
    

}

chr_t::~chr_t()
{

}

} /* namespace ccat_aux */
