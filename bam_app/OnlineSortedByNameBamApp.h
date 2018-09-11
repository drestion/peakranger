/*
 * OnlineSortedByNameBamApp.h
 *
 *  Created on: May 29, 2012
 *      Author: xfeng
 */

#ifndef ONLINESORTEDBYNAMEBAMAPP_H_
#define ONLINESORTEDBYNAMEBAMAPP_H_
#include <string>
#include <ostream>
#include <stdint.h>
#include "bam_app/OnlineBamAppImp.h"
namespace bam_app {

/*
 * A framework for apps that scale with
 * two adjacent bam reads with the
 * same read name, using a single
 * bam file. The bam file is assumed
 * to be sorted by read name.
 * samtools sort -n sorts the file.
 */

class OnlineSortedByNameBamApp {
public:
    OnlineSortedByNameBamApp();
    virtual ~OnlineSortedByNameBamApp();
    OnlineSortedByNameBamApp(aux::OnlineBamAppImp* parser);

    virtual void processReads(const std::string& file, std::ostream& os);
    virtual void report(std::ostream& os);
};

} /* namespace bam_app */
#endif /* ONLINESORTEDBYNAMEBAMAPP_H_ */
