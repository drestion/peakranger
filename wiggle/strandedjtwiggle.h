/*
 * strandedjtwiggle.h
 *
 *  Created on: Jan 11, 2012
 *      Author: xfeng
 */

#ifndef STRANDEDJTWIGGLE_H_
#define STRANDEDJTWIGGLE_H_

#include "JTwigglefile.h"

class stranded_jtwiggle: public JT_wiggle_file {
public:
    stranded_jtwiggle();
    virtual ~stranded_jtwiggle();
    virtual void export_wiggle(Reads& reads, std::ostream& os);

    virtual void export_wiggle(std::vector<uint32_t>& preads,
            std::vector<uint32_t>& nreads, std::string chr, std::ostream& os);
    virtual void export_wiggle_pos(std::vector<uint32_t>& preads, std::string chr,
            std::ostream& os);
    virtual void export_wiggle_neg(std::vector<uint32_t>& nreads, std::string chr,
            std::ostream& os);
    virtual void export_wiggle(Reads& reads, const char* file);
    virtual void split_export_wiggle(Reads& reads, std::ostream& os);
    virtual void split_export_wiggle(Reads& reads, const char* file);
    virtual void export_wiggle_gzip(Reads& reads, const char* file);

    virtual void split_export_wiggle_gzip(Reads& reads, const char* file);
protected:
    virtual void _process_neg(uint32_t start, uint32_t end, uint32_t readlength,
            uint32_t readextlength, std::vector<uint32_t>::iterator nreadsStart,
            std::vector<uint32_t>::iterator nreadsEnd, std::ostream& os);

    std::vector<uint32_t> _ncolorRGB;
};

#endif /* STRANDEDJTWIGGLE_H_ */
