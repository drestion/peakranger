/*
 * peakranger_cmd_option_parser
 *
 *  Created on: July 25, 2011
 *      Author: xin
 */

#ifndef PEAKRANGEROPTION_H_
#define PEAKRANGEROPTION_H_

#include <iostream>
#include "option_parser/cmd_option_parser.h"
#include "option_parser/OptionAux.h"

class peakranger_cmd_option_parser: public cmd_option_parser {
public:
    peakranger_cmd_option_parser(int argc, char** argv,const std::string& version);

    ~peakranger_cmd_option_parser();

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

#endif /* PEAKRANGEROPTION_H_ */
