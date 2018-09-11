/*
 * LC.h
 *
 *  Created on: Jul 10, 2012
 *      Author: xfeng
 */

#ifndef LC_H_
#define LC_H_
#include <string>
namespace app {

class LC {
public:
    LC();
    virtual ~LC();
    static void run(int argc, char** argv);
    static std::string version;
};

} /* namespace app */
#endif /* LC_H_ */
