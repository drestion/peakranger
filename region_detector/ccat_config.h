/*
 * ccat_config.h
 *
 *  Created on: Apr 3, 2012
 *      Author: xin
 */

#ifndef CCAT_CONFIG_H_
#define CCAT_CONFIG_H_

#include <string>
#include <vector>
#include <iostream>
namespace ccat_aux {

class ccat_config_t {
    friend std::ostream& operator<<(std::ostream& os,
            const ccat_config_t& conf) {
        os << "fragmentSize:" << conf.fragmentSize;
        os << "\nslidingWinSize:" << conf.slidingWinSize;
        os << "\nmovingStep:" << conf.movingStep;
        os << "\nminCount:" << conf.minCount;
        os << "\nisStrandSensitiveMode:" << conf.isStrandSensitiveMode;
        os << "\noutputNum:" << conf.outputNum;
        os << "\nrandomSeed:" << conf.randomSeed;
        os << "\nminScore:" << conf.minScore;
        os << "\nbootstrapPass:" << conf.bootstrapPass;
        os << "\nsmoothingFactor:" << conf.smoothingFactor;
        return os;
    }
public:
    ccat_config_t();
    virtual ~ccat_config_t();
    int fragmentSize;
    int slidingWinSize;
    int movingStep;
    int minCount;
    int isStrandSensitiveMode;
    int outputNum;
    int randomSeed;
    double minScore;
    int bootstrapPass;
    double smoothingFactor;

};

} /* namespace ccat_aux */
#endif /* CCAT_CONFIG_H_ */
