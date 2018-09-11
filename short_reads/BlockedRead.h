/*
 * BlockedRead.h
 *
 *  Created on: May 10, 2012
 *      Author: xin
 */

#ifndef BLOCKEDREAD_H_
#define BLOCKEDREAD_H_

#include "common/stl_header.h"
#include "common/boost_header.h"
#include "Read.h"

namespace reads {
/*
 * Assume all blocks are on the same chromosome
 */
class BlockedRead {
    friend std::ostream& operator<<(std::ostream& os, const BlockedRead& rhs) {
        foreach(Read r, rhs.mBlocks) {
            os << r;
        }
        return os;
    }

public:
    BlockedRead();
    virtual ~BlockedRead();

    BlockedRead(const std::vector<Read>& blocks);

    void getBlocks(std::vector<Read>& res) const;
    std::vector<Read> getBlocks() const;
    std::string getChr() const;

    void setBlocks(const std::vector<Read>& blocks) {
        mBlocks = blocks;
    }
    Strand getDir() const;

    int32_t getStart() const;

    bool operator==(const BlockedRead& rhs) const {
        if (this->mBlocks.size() != rhs.mBlocks.size()) {
            return false;
        } else {
            return std::equal(mBlocks.begin(), mBlocks.end(),
                    rhs.mBlocks.begin());
        }
    }
protected:
    std::vector<Read> mBlocks;

};

void offset(BlockedRead& read, int32_t length);
}
#endif /* BLOCKEDREAD_H_ */
