/*
 * WigPE.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef WIGPE_H_
#define WIGPE_H_
#include <string>
namespace app {

class WigPE {
public:
    WigPE();
    virtual ~WigPE();
    static void run(int argc, char** argv);
    static std::string version;

private:
    static std::string normOutFileName(const std::string& fn);
};

} /* namespace app */
#endif /* WIGPE_H_ */
