/*
 * stringutil.h
 *
 *  Created on: Oct 6, 2011
 *      Author: xfeng
 */

#ifndef STRINGUTIL_H_
#define STRINGUTIL_H_

#include <string>
#include <boost/algorithm/string.hpp>
#include <sstream>
namespace utils {

class stringutil {
public:
    stringutil();
    virtual ~stringutil();

    static void remove_last(std::string& s, const char* toremove);
    static void get_dir_file(std::string s, std::string& dir, std::string& file,
            std::string& ext);
};

template<typename T>
std::string vector_to_string(const T& vec, const char* sep) {
    std::stringstream str;

    typename T::const_iterator it;
    it = vec.begin();
    for (; it != vec.end(); it++) {
        str << *it << sep;
    }
    std::string result(str.str());
    boost::replace_last(result, sep, "");
    return result;
}

template<typename T>
std::string vector_to_string(const T& vec, std::string sep) {
    std::stringstream str;

    typename T::const_iterator it;
    it = vec.begin();
    for (; it != vec.end(); it++) {
        str << *it << sep;
    }
    std::string result(str.str());
    boost::replace_last(result, sep, "");
    return result;
}

/*
 * based on: http://www.cplusplus.com/forum/beginner/1962/
 */
std::string ExtractDirectory(const std::string& path);
/*
 * based on: http://www.cplusplus.com/forum/beginner/1962/
 */
std::string ExtractFilename(const std::string& path);
std::string ExtractExt(const std::string& path);
/*
 * based on: http://www.cplusplus.com/forum/beginner/1962/
 */
std::string ChangeExtension(const std::string& path, const std::string& ext);
}
#endif /* STRINGUTIL_H_ */
