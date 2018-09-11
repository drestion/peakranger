/*
 * Runnable.h
 *
 *  Created on: Jul 18, 2012
 *      Author: xfeng
 */

#ifndef RUNNABLE_H_
#define RUNNABLE_H_

namespace threads {
/**
 * Implement this class to use
 * threads facility such as
 * IndependentRunner
 */
template<typename Job>
class Runnable {
public:
    Runnable(){};
    virtual ~Runnable(){};

    /**
     * The details on running the job
     */
    virtual void runJob(Job& job) = 0;
};

}


 /* namespace threads */
#endif /* RUNNABLE_H_ */
