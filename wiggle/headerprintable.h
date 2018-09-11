/*
 * headerprintable.h
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#ifndef HEADERPRINTABLE_H_
#define HEADERPRINTABLE_H_
#include <ostream>
#include <vector>
#include <string>
#include <stdint.h>

class header_printable {
public:
    header_printable();
    virtual ~header_printable();
    void print_wigfile_trackheader(std::ostream& pof, std::vector<uint32_t> col,
            const char* _name = "", uint32_t _priority = 20) {
        if (col.size() < 3) {
            uint32_t ga[3] = { 50, 126, 184 };
            col = std::vector<uint32_t>(ga, ga + 3);
        }
        pof << "track type=wiggle_0 name=\"" << _name << "\" "
                << "visibility=dense " << "color=" << col[0] << "," << col[1]
                << "," << col[2] << " " << "altColor=" << col[0] << ","
                << col[1] << "," << col[2] << " " << "priority=" << _priority
                << "\n";
    }
};

#endif /* HEADERPRINTABLE_H_ */
