/*
 * util_print.h
 *
 *  Created on: Jan 5, 2012
 *      Author: xfeng
 */

#ifndef UTIL_PRINT_H_
#define UTIL_PRINT_H_

#include <string>
#include <iostream>

namespace utilprint {
class Msg {
public:
    virtual ~Msg(){}
    virtual std::string tostring() {
        return _msg;
    }
    void set_msg(const std::string& msg) {
        _msg = msg;
    }
    void set_msg(const char* msg) {
        _msg = std::string(msg);
    }
    void print_msg(std::ostream& os,const char* del="#") {
        os << del <<_msg;
    }
protected:
    std::string _msg;
};

class citation: public Msg {
public:
    citation() {
        _msg.append("Feng X, Grossman R, Stein L: PeakRanger:");
        _msg.append("A cloud-enabled peak caller for ChIP-seq data.");
        _msg.append("BMC Bioinformatics 2011, 12(1):139.\n");
    }

};

class ccatCitation: public Msg {
public:
    ccatCitation() {
        _msg.append("Xu, H., L. Handoko, et al. (2010).");
        _msg.append("A signal-noise model for significance ");
        _msg.append("analysis of ChIP-seq with negative control.");
        _msg.append("Bioinformatics 26(9): 1199-1204.\n");
    }

};

class bcpCitation: public Msg {
public:
    bcpCitation() {
    	_msg.append("Xing H, Mo Y, Liao W, Zhang MQ (2012) Genome-Wide Localization of");
    	_msg.append(" Protein-DNA Binding and Histone Modification by a Bayesian");
    	_msg.append(" Change-Point Method with ChIP-seq Data. PLoS Comput Biol 8(7): e1002613.");
    	_msg.append(" doi:10.1371/journal.pcbi.1002613\n");

    }

};
const std::string HEADER("\033[95m");
const std::string OKBLUE("\033[94m");
const std::string OKGREEN("\033[92m");
const std::string WARNING("\033[93m");
const std::string SECTION("\033[93m");
const std::string FAIL("\033[91m");
const std::string ENDC("\033[0m");
}
#endif /* UTIL_PRINT_H_ */
