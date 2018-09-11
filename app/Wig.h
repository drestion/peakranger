/*
 * Wig.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef WIGAPP_H_
#define WIGAPP_H_
#include <string>
namespace app {

class Wig {
public:
    Wig();
    virtual ~Wig();
    static void run(int argc, char** argv);
    static std::string version;
};

} /* namespace app */
#endif /* WIGAPP_H_ */
