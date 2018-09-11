/*
 * CCAT.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef CCATAPP_H_
#define CCATAPP_H_
#include <string>
namespace app {

class CCAT {
public:
    CCAT();
    virtual ~CCAT();
    static void run(int argc, char** argv);
    static std::string version;
};

} /* namespace app */
#endif /* CCATAPP_H_ */

