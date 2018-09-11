/*
 * exceptions.h
 *
 *  Created on: Apr 28, 2011
 *      Author: xin
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_
#include <string>
#include <stdio.h>
#include <stdexcept>

class RangerException{
public:
    virtual ~RangerException() {
    }
    RangerException(std::string& msg) :
            _msg(msg) {
    }

    RangerException(const char* msg) :
            _msg(msg) {
    }
    RangerException() :
            _msg() {
    }
    ;
    std::string error() {
        return _msg;
    }
    virtual void debugPrint() {
        printf("default exception msg:%s", _msg.c_str());
    }

protected:
    std::string _msg;
};
struct QuerySeqIsEmpty: public RangerException {
    QuerySeqIsEmpty(const char* msg) :
            RangerException(msg) {
    }
};
struct IllegalCigarString: public RangerException {
    IllegalCigarString(const char* msg) :
            RangerException(msg) {
    }
};
class NotPairEndBamRead: public RangerException {
public:
    NotPairEndBamRead(const char* msg) :
            RangerException(msg) {
    }
};
class FileNotReady: public RangerException {
public:

    FileNotReady(const char* msg) :
            RangerException(msg) {
    }
    FileNotReady(char* msg) :
            RangerException(msg) {
    }
    FileNotReady(std::string& msg) :
            RangerException(msg) {
    }
    void debugPrint() {
        printf("File not ready:\n%s", _msg.c_str());
    }
};

class FileNotGood: public RangerException {
public:
    FileNotGood() {
    }
    FileNotGood(const char* msg) :
            RangerException(msg) {

    }

    FileNotGood(std::string& msg) :
            RangerException(msg) {
    }
    virtual void debugPrint() {
        printf("Can not open this file: %s\n", _msg.c_str());
    }
};

class ChrNotFound: public RangerException {
public:
    ChrNotFound(std::string& chr) :
            RangerException(chr) {
    }
    ChrNotFound(const char* chr) :
            RangerException(chr) {
    }

    void debugPrint() {
        printf("The chromosome: %s is not valid.", _msg.c_str());
    }
};

class DataLineNotValid: public RangerException {
public:
    DataLineNotValid(std::string& line, std::string& expectedFormat) :
            RangerException(line), _format(expectedFormat) {
    }
    DataLineNotValid(std::string& line) :
            RangerException(line), _format("NOT SPECIFIED") {
    }
    void debugPrint() {
        printf("Found an invalid line in the data file that violates "
                "the specified %s format:\n %s ", _format.c_str(),
                _msg.c_str());
    }
private:
    std::string _format;
};

class InvalidArgument: public RangerException {
public:
    InvalidArgument(const char* msg) :
            RangerException(msg) {
    }
};

class DebugFlag: public RangerException {
public:
    DebugFlag(const char* msg) :
            RangerException(msg) {
    }
    DebugFlag(std::string& msg) :
            RangerException(msg) {
    }
};

class not_in_range: public RangerException {
public:
    not_in_range(const char* msg) :
            RangerException(msg) {
    }
};

#endif /* EXCEPTIONS_H_ */
