/*
 * ImgScriptRunner.cpp
 *
 *  Created on: Jul 18, 2012
 *      Author: xfeng
 */

#include "ggplay/ImgScriptRunner.h"
#include "ggplay/LinuxRCmd.h"
#include "threading/IndependentRunner.h"

using namespace threads;

namespace ggplay {

ImgScriptRunner::ImgScriptRunner(const std::vector<std::string>& mScripts,
        uint32_t threads) :
        mScripts(mScripts), mThreads(threads) {
}

ImgScriptRunner::~ImgScriptRunner() {
}

void ImgScriptRunner::run() {
    LinuxRCmd r;
    IndependentRunner<std::string, LinuxRCmd> scriptCPU(mScripts, r, mThreads);
    scriptCPU.run();
}

} /* namespace ggplay */
