/*
 * WigPEOption.h
 *
 *  Created on: Jun 25, 2012
 *      Author: xfeng
 */

#ifndef WIGPEOPTION_H_
#define WIGPEOPTION_H_
#include <boost/program_options.hpp>
#include "option_parser/OptionAux.h"
namespace options {

namespace po = boost::program_options;

class WigPEOption {
public:
    WigPEOption(const std::string& version);
    virtual ~WigPEOption();
    void parse(int argc, char** argv);

    std::string printParsedOpts();
    std::string version;
    static int min_args;
    boost::program_options::variables_map mVM;
private:
    void hasEnoughArgs(int argc);
    void verifyOptions();
    po::opt all;
    po::p_opt popt;
};

} /* namespace options */
#endif /* WIGPEOPTION_H_ */
