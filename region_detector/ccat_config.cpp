/*
 * ccat_config.cpp
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#include "ccat_config.h"

namespace ccat_aux {

ccat_config_t::ccat_config_t()
: fragmentSize(0),
  slidingWinSize(0),
  movingStep(0),
  minCount(0),
  isStrandSensitiveMode(0),
  outputNum(0),
  randomSeed(0),
  minScore(0),
  bootstrapPass(0),
  smoothingFactor(0)
{

}

ccat_config_t::~ccat_config_t()
{

}

} /* namespace ccat_aux */
