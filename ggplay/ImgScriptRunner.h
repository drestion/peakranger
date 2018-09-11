/*
 * ImgScriptRunner.h
 *
 *  Created on: Jul 18, 2012
 *      Author: xfeng
 */

#ifndef IMGSCRIPTRUNNER_H_
#define IMGSCRIPTRUNNER_H_
#include <string>
#include <vector>
namespace ggplay {

class ImgScriptRunner {
public:
    ImgScriptRunner(const std::vector<std::string>& mScripts, uint32_t threads);
    virtual ~ImgScriptRunner();

    void run();

private:
    std::vector<std::string> mScripts;
    uint32_t mThreads;
};

} /* namespace ggplay */
#endif /* IMGSCRIPTRUNNER_H_ */
