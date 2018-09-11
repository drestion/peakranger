/*
 * BCPOption.h
 *
 *  Created on: Sep 27, 2014
 *      Author: Xin Feng
 */

#ifndef BCPOPTION_H_
#define BCPOPTION_H_
#include <boost/program_options.hpp>
#include "option_parser/cmd_option_parser.h"

namespace options {
namespace po = boost::program_options;
class BCPOption : public cmd_option_parser{
public:
    BCPOption();
    virtual ~BCPOption();
    void parse();
    BCPOption(int argc, char** argv, const std::string& version);
    std::string printParsedOpts() const;

    static int min_args;
    std::string version;
    po::variables_map mVM;
    bool needHtml() const;
    void setNeedHtml(bool _html);
    void print_option(std::ostream& os);
       void print_option_file(std::ostream& os) const;

private:
    void hasEnoughArgs(int argc);
    void verifyOptions();
    bool _html;
    po::opt all;
    po::p_opt popt;
    po::opt running_modes;
};

} /* namespace options */
#endif /* BCPOPTION_H_ */
