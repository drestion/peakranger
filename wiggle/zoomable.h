/*
 * zoomable.h
 *
 *  Created on: Jan 25, 2012
 *      Author: xfeng
 */

#ifndef ZOOMABLE_H_
#define ZOOMABLE_H_

#include "wig.h"
#include "wigbuilder.h"
#include <stdint.h>

class zoomable{
public:
    zoomable();
    virtual ~zoomable();

    virtual void smooth(wigs& _wigs,
                        wigs& result,
                        uint32_t window,
                        uint32_t overlap) = 0;
};

#endif /* ZOOMABLE_H_ */
