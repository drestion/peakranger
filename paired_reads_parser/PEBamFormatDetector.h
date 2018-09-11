/*
 * PEBamFormatDetector.h
 *
 *  Created on: Jun 6, 2012
 *      Author: xfeng
 */

#ifndef PEBAMFORMATDETECTOR_H_
#define PEBAMFORMATDETECTOR_H_

#include <string>
#include <stdint.h>
namespace parser {

class PEBamFormatDetector {
public:
    PEBamFormatDetector();
    virtual ~PEBamFormatDetector();
    virtual bool isPairEndBamFile(const std::string& file,uint32_t linesToTest);
};

} /* namespace parser */
#endif /* PEBAMFORMATDETECTOR_H_ */
