/*
 * chip4cwigcmdoption.h
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#ifndef CHIP4CWIGCMDOPTION_H_
#define CHIP4CWIGCMDOPTION_H_

#include "wig_cmd_option.h"
#include <iostream>
class chip4c_wig_cmd_option: public wig_cmd_option {
public:
    chip4c_wig_cmd_option(int argc, char** argv) :
            wig_cmd_option(argc, argv,"1.16") {

    }

    void parse();
    void print_option(std::ostream& os);
    void printHelp() const;
    bool isSplit() const;
    void setSplit(bool _split);
    bool isGz() const;
    bool isStranded() const;
    void setGz(bool _gz);
    void setStranded(bool _stranded);
    void print_option_file(std::ostream & os) const;
    uint32_t getOverlapSz() const;
    uint32_t getWinwdowSz() const;
    void setOverlapSz(uint32_t overlapSz);
    void setWinwdowSz(uint32_t winwdowSz);
    uint32_t getRepeat() const;
    void setRepeat(uint32_t _repeat);
    std::string version;
private:

    po::variables_map vm;
    void verify();
    inline void require(const char* opt, po::variables_map& vm);
    bool _split;
    bool _gz;
    bool _stranded;
    uint32_t _window_sz;
    uint32_t _overlap_sz;
    uint32_t _repeat;
};

#endif /* CHIP4CWIGCMDOPTION_H_ */
