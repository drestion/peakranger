/*
 * ranger4ccmdparser.h
 *
 *  Created on: Feb 7, 2012
 *      Author: xfeng
 */

#ifndef RANGER4CCMDPARSER_H_
#define RANGER4CCMDPARSER_H_

#include "cmd_option_parser.h"
#include "utils/exceptions.h"
#include <stdlib.h>
#include <boost/program_options.hpp>
#include <iostream>

class ranger4c_cmd_parser: public cmd_option_parser {
public:
    ranger4c_cmd_parser(int argc, char** argv) :
            cmd_option_parser(argc, argv) {
    }

    virtual ~ranger4c_cmd_parser() {
    }

    void parse();
    void print_option(std::ostream& os);
    void print_option_file(std::ostream& os) const;
    void printHelp() const;
    bool isSplit() const;
    void setSplit(bool _split);
    uint32_t getW() const;
    void setW(uint32_t _w);
    bool isHas4CControl() const;
    void setHas4CControl(bool _has4CControl);
    std::string version;
private:
    bool _split;
    bool _has4CControl;
    uint32_t _w;

    po::variables_map vm;

    void verify();

};

#endif /* PEAKRANGEROPTION_H_ */
