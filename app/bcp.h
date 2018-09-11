/*
 * bcp.h
 *
 *  Created on: Sep 27, 2014
 *      Author: Xin Feng
 */

#ifndef BCP_H_
#define BCP_H_ 
#include <string>
namespace app {

class BCP {
public:
    BCP();
    virtual ~BCP();
    static void run(int argc, char** argv);
    static std::string version;
};

} /* namespace app */
#endif /* BCP_H_ */

