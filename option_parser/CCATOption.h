/*
 * CCATOption.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef CCATOPTION_H_
#define CCATOPTION_H_
#include <boost/program_options.hpp>
#include "option_parser/OptionAux.h"

namespace options {
namespace po = boost::program_options;
class CCATOption {
public:
    CCATOption();
    virtual ~CCATOption();
    void parse(int argc, char** argv);

    std::string printParsedOpts();

    static int min_args;
    std::string version;
    boost::program_options::variables_map mVM;
private:
    void hasEnoughArgs(int argc);
    void verifyOptions();
    po::opt all;
    po::p_opt popt;
    po::opt running_modes;
};

} /* namespace options */
#endif /* CCATOPTION_H_ */
