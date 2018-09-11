/*
 * JTwigglefile.h
 *
 *  Created on: Oct 4, 2011
 *      Author: xfeng
 */

#ifndef JTWIGGLEFILE_H_
#define JTWIGGLEFILE_H_

#include "wiggle_reporter.h"

/*
 * The difference between this one and its mother is the
 * way the neg reads were extended.
 *
 * chr1 start end -
 *
 * is extended to:
 *
 * start = start + readlength -extlength
 * end = end = start + readlength;
 *
 * the way it works is by reimplementing _process
 *
 */
class JT_wiggle_file: public wiggle_reporter {
public:
    JT_wiggle_file();
    virtual ~JT_wiggle_file();
    virtual void export_wiggle(Reads& reads, std::ostream& os);
    virtual void export_wiggle(std::vector<uint32_t>& preads,
            std::vector<uint32_t>& nreads, std::string chr, std::ostream& os);
    virtual void export_wiggle(Reads& reads, const char* file);
    virtual void split_export_wiggle(Reads& reads, std::ostream& os);
    virtual void split_export_wiggle(Reads& reads, const char* file);
    virtual void export_wiggle_gzip(Reads& reads, const char* file);

    virtual void split_export_wiggle_gzip(Reads& reads, const char* file);
protected:
    virtual void _process(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator readsStart,
            std::vector<uint32_t>::iterator readsEnd,
            std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, std::ostream& os);
};

#endif /* JTWIGGLEFILE_H_ */
