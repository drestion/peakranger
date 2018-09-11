/*
 * LinuxRCmd.h
 *
 *  Created on: Jul 18, 2012
 *      Author: xfeng
 */

#ifndef LINUXRCMD_H_
#define LINUXRCMD_H_

#include "threading/Runnable.h"
#include <string>
namespace ggplay {

class LinuxRCmd: public threads::Runnable<std::string> {
public:
    LinuxRCmd();
    virtual ~LinuxRCmd();

    void runJob(std::string& RCmd);

private:
    int rmScript(std::string& rf);
};

} /* namespace ggplay */
#endif /* LINUXRCMD_H_ */
