/*
 * PeakRanger.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef PeakRanger_H_
#define PeakRanger_H_
#include <string>
namespace app {

class PeakRanger {
public:
    PeakRanger();
    virtual ~PeakRanger();
    static void run(int argc, char** argv);
    static std::string version;
};

} /* namespace app */
#endif /* PeakRanger_H_ */
