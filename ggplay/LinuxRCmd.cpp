/*
 * LinuxRCmd.cpp
 *
 *  Created on: Jul 18, 2012
 *      Author: xfeng
 */

#include <stdlib.h>
#include "ggplay/LinuxRCmd.h"
using namespace std;
namespace ggplay {

LinuxRCmd::LinuxRCmd() {

}

LinuxRCmd::~LinuxRCmd() {

}

int LinuxRCmd::rmScript(std::string& rf) {
    string rm = "rm -f " + rf;
    return system(rm.c_str());
}

void LinuxRCmd::runJob(std::string& rf) {
    string R = "R CMD BATCH " + rf + " /dev/null";
    int r1 = system(R.c_str());
    r1 += rmScript(rf);
}

} /* namespace ggplay */
