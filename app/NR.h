/*
 * NR.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef NR_H_
#define NR_H_
#include <string>
namespace app {

class NR {
public:
    NR();
    virtual ~NR();
    static void run(int argc, char** argv);
    static std::string version;
};

} /* namespace app */
#endif /* NR_H_ */
