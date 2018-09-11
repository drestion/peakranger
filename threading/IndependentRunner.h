/*
 * IndependentRunner.h
 *
 *  Created on: Jul 18, 2012
 *      Author: xfeng
 */

#ifndef INDEPENDENTRUNNER_H_
#define INDEPENDENTRUNNER_H_
#include <boost/thread.hpp>
#include <boost/bind.hpp>

namespace threads {
/**
 * Runs a vector of jobs that can run without
 * depending on other jobs
 */
template<typename Job, typename JobRunner>
class IndependentRunner {
public:
    IndependentRunner(std::vector<Job>& jobs, JobRunner& runner,
            uint32_t threads);
    ~IndependentRunner();

public:
    void run();

protected:
    void runJob(Job& job);
    void runJobs();

private:
    boost::mutex mMutex;
    uint32_t mThreads;
    std::vector<Job>& mJobs;
    JobRunner& mRunner;
};

template<typename Job, typename JobRunner>
IndependentRunner<Job, JobRunner>::IndependentRunner(std::vector<Job>& jobs,
        JobRunner& runner, uint32_t threads) :
        mMutex(), mThreads(threads), mJobs(jobs), mRunner(runner) {
}

template<typename Job, typename JobRunner>
IndependentRunner<Job, JobRunner>::~IndependentRunner() {
}

template<typename Job, typename JobRunner>
void IndependentRunner<Job, JobRunner>::run() {
 //   boost::thread_group threads;
//    size_t i = 0;
//    while (i++ < mThreads) {
  //     threads.create_thread(boost::bind(&IndependentRunner::runJobs, this));
  //  }
   // threads.join_all();
   runJobs();
}

template<typename Job, typename JobRunner>
void IndependentRunner<Job, JobRunner>::runJob(Job& job) {
    mRunner.runJob(job);
}

template<typename Job, typename JobRunner>
void IndependentRunner<Job, JobRunner>::runJobs() {
    Job job;
    while (true) {
        {
            boost::mutex::scoped_lock l(mMutex);
            if (mJobs.empty()) {
                break; //the exit
            }
            job = mJobs.back();
            mJobs.pop_back();
        }
        runJob(job);
    }
}

}
/* namespace threads */
#endif /* INDEPENDENTRUNNER_H_ */
