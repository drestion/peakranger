/*
 * unequalstrandwigglereporter.h
 *
 *  Created on: Sep 27, 2011
 *      Author: xfeng
 */

#ifndef UNEQUALSTRANDWIGGLEREPORTER_H_
#define UNEQUALSTRANDWIGGLEREPORTER_H_

#include "wiggle_reporter.h"

/*
 * This class functions similarily with wiggle_reporter.h
 * However, it can handle the case when the reads are missing
 * in one strand of a chromosome
 */
class unequal_strand_wiggle_reporter:public wiggle_reporter{
public:
    unequal_strand_wiggle_reporter();
    virtual ~unequal_strand_wiggle_reporter();

    virtual void split_export_wiggle(Reads& reads,
                                     const char* file);
    virtual void export_wiggle(Reads& reads,
                               std::ostream& os);
    virtual void export_wiggle(Reads& reads,
                               const char* file);
};

#endif /* UNEQUALSTRANDWIGGLEREPORTER_H_ */
