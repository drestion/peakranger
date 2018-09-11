/*
 * CigarStringFactory.h
 *
 *  Created on: May 12, 2012
 *      Author: xfeng
 */

#ifndef CIGARSTRINGFACTORY_H_
#define CIGARSTRINGFACTORY_H_

#include <boost/shared_ptr.hpp>
#include "bamtools/BamAux.h"
#include "CigarString.h"
#include <stdint.h>
namespace reads {

class CigarStringFactory {
public:
	static boost::shared_ptr<CigarString> buildCigar(
			const BamTools::CigarOp& cigar, const int32_t offset);

	static boost::shared_ptr<CigarString> buildCigar(
			const char& ch, const int32_t offset, const int32_t length);
};

} /* namespace reads */
#endif /* CIGARSTRINGFACTORY_H_ */
