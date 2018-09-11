/*
 * ccatcmdoptionparser.h
 *
 *  Created on: Jun 28, 2012
 *      Author: xfeng
 */

#ifndef CCATCMDOPTIONPARSER_H_
#define CCATCMDOPTIONPARSER_H_
#include "option_parser/cmd_option_parser.h"
#include "option_parser/OptionAux.h"
#include <iostream>
namespace options {

class ccat_cmd_option_parser: public cmd_option_parser {
public:
    ccat_cmd_option_parser(int argc, char** argv, const std::string& version);

    ~ccat_cmd_option_parser();

    void parse();
    void print_option(std::ostream& os);
    void print_option_file(std::ostream& os) const;

    bool isSplit() const;
    void setSplit(bool _split);
    bool needHtml() const;
    void setNeedHtml(bool _html);
    std::string version;

private:
    bool _split;
    bool _html;
    po::variables_map vm;
    po::opt all;
    po::p_opt popt;

    void verify();

};

} /* namespace options */
#endif /* CCATCMDOPTIONPARSER_H_ */
