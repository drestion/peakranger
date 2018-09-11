/*
 * OptionAux.h
 *
 *  Created on: Jun 25, 2012
 *      Author: xfeng
 */

#ifndef OPTIONAUX_H_
#define OPTIONAUX_H_
#include <boost/program_options.hpp>
#include <iostream>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;
namespace boost {
namespace program_options {
typedef options_description opt;
typedef positional_options_description p_opt;
}
}

namespace options {
namespace aux {

class SharedOpts {

    SharedOpts();
    static boost::program_options::opt help_verbose_version;
};
void printHelp(const boost::program_options::options_description& opts,
        std::ostream& os);
void file_r_good(const char* file);
void file_w_good(const char* file);
void printVersion(std::ostream& os, const std::string& version);
void require(const char* opt, boost::program_options::variables_map& vm);
void requireNotEqual(const char* opt, const char* notequal, const char* exepctedvalue,
        boost::program_options::variables_map& vm);
template<typename T>
bool is_in_range(T val, T low, T high) {
    if (val > low && val < high)
        return true;
    return false;
}

template<typename T>
std::string vector_to_string(std::vector<T>& vec, uint32_t width_count,
        std::string prefix) {
    std::stringstream str;
    typename std::vector<T>::iterator it;
    it = vec.begin();
    uint32_t cnt = 0;
    for (; it != vec.end(); it++) {
        str << *it << ",";
        cnt++;
        if (cnt == width_count) {
            str << "\n";
            size_t prefixlength = prefix.size();
            while (prefixlength-- > 0) {
                str << " ";
            }
            cnt = 0;
        }
    }
    std::string result(str.str());
    boost::replace_last(result, ",", "");
    return result;
}
void parseChrTable(std::string& filename, std::vector<std::string>& results);
std::string printRawOpts(int argc, char** argv);
}
} /* namespace options */
#endif /* OPTIONAUX_H_ */
