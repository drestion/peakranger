/*
 * Guarded.h
 *
 *  Created on: May 13, 2012
 *      Author: xfeng
 */

#ifndef GUARDED_H_
#define GUARDED_H_

#include <string>

namespace utils {

template<typename E>
class Guarded {
public:

    Guarded(bool cond, const char* msg) :
            mMsg(std::string(msg)) {
        if (cond) {
            throw E(msg);
        }
    }
    virtual ~Guarded() {
    }
    ;

private:
    std::string mMsg;
};

} /* namespace utils */
#endif /* GUARDED_H_ */
