/*
 * NROption.h
 *
 *  Created on: Jun 27, 2012
 *      Author: xfeng
 */

#ifndef NROPTION_H_
#define NROPTION_H_
#include <boost/program_options.hpp>
#include "option_parser/OptionAux.h"
namespace options {
namespace po = boost::program_options;
class NROption {
public:
    NROption(const std::string& version);
    virtual ~NROption();
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
#endif /* NROPTION_H_ */
