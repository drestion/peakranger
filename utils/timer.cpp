/*
 * timer.cpp
 *
 *  Created on: Jun 28, 2012
 *      Author: xfeng
 */

#include "utils/timer.h"
using namespace std;
namespace utils {
void logDate(std::ostream& os, bool nl) {
    time_t t = time(0); // get time now
    struct tm * now = localtime(&t);
    os << "#";
    os << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-'
            << now->tm_mday << " " << now->tm_hour - 1 << ":" << now->tm_min
            << ":" << now->tm_sec;
    if (nl)
        os << std::endl;
}

void getDate(std::string& date) {
    stringstream os;
    time_t t = time(0); // get time now
    struct tm * now = localtime(&t);

    os << (now->tm_year + 1900) << '-' << (now->tm_mon + 1) << '-'
            << now->tm_mday;

    date = os.str();
}
void logTime(std::ostream& os, bool nl) {
    struct tm *current;
    time_t now;
    time(&now);
    current = localtime(&now);
    os << setfill('0') << setw(2) << current->tm_hour << ":" << setfill('0')
            << setw(2) << current->tm_min << ":" << setfill('0') << setw(2)
            << current->tm_sec;
//    os <<"\n";
    if (nl)
        os << std::endl;
}
}
