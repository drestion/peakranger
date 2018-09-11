/*
 * stringutil.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: xfeng
 */

#include "stringutil.h"
#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
using namespace std;

namespace utils {

stringutil::stringutil() {

}

stringutil::~stringutil() {

}

void stringutil::remove_last(string & s, const char* toremove) {

    size_t ind = s.find_last_of(toremove);

    if (ind != string::npos) {
        s = s.substr(0, ind);
    }
}

void stringutil::get_dir_file(string s, string& dir, string& file,
        string& ext) {
    string str = s;
    size_t found;
    found = str.find_last_of("/\\");
    if (found != string::npos) {
        dir = str.substr(0, found);
        file = str.substr(found + 1);
    } else {
        dir = ".";
        file = str;
    }
    found = str.find_last_of(".");
    if (found != string::npos) {
        if (boost::starts_with(str, "./") || boost::starts_with(str, "")) {
            string afters = str.substr(str.find_first_of("/") + 1);
            found = afters.find_last_of(".");
            ext = afters.substr(found + 1);
        } else {
            ext = str.substr(found + 1);
        }
    } else {
        ext = "";
    }
}

std::string ExtractDirectory(const std::string& path) {
    return path.substr(0, path.find_last_of("/\\") + 1);
}

std::string ExtractFilename(const std::string& path) {
    return path.substr(path.find_last_of("/\\") + 1);
}

std::string ExtractExt(const std::string& path) {
    std::string filename = ExtractFilename(path);
    return filename.substr(filename.find_last_of('.') + 1);
}

std::string ChangeExtension(const std::string& path, const std::string& ext) {
    std::string filename = ExtractFilename(path);
    return ExtractDirectory(path)
            + filename.substr(0, filename.find_last_of('.')) + ext;
}
}
