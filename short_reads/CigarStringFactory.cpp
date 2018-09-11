/*
 * CigarStringFactory.cpp
 *
 *  Created on: May 12, 2012
 *      Author: xfeng
 */

#include "CigarStringFactory.h"
#include "common/ranger_debug.h"
#include <boost/make_shared.hpp>
using namespace boost;

namespace reads {

enum {
	cigarM = 'M',
	cigarI = 'I',
	cigarD = 'D',
	cigarN = 'N',
	cigarS = 'S',
	cigarH = 'H',
	cigarP = 'P'
};

boost::shared_ptr<CigarString> reads::CigarStringFactory::buildCigar(
		const BamTools::CigarOp& c, const int32_t offset) {
	return buildCigar(c.Type, offset, c.Length);
}

boost::shared_ptr<CigarString> CigarStringFactory::buildCigar(const char& c,
		const int32_t offset,const int32_t length) {
	shared_ptr<CigarString> res;
	switch (c) {
	case cigarM:
		res = make_shared<CigarM>(length, offset);
		break;
	case cigarI:
		res = make_shared<CigarI>(length, offset);
		break;
	case cigarD:
		res = make_shared<CigarD>(length, offset);
		break;
	case cigarN:
		res = make_shared<CigarN>(length, offset);
		break;
	case cigarS:
		res = make_shared<CigarS>(length, offset);
		break;
	case cigarH:
		res = make_shared<CigarH>(length, offset);
		break;
	case cigarP:
		res = make_shared<CigarP>(length, offset);
		break;
	default:
		std::string str("Illegal cigar type: ");
		str.append(1, c);
		throw IllegalCigarString(str.c_str());
	}
	return res;
}

}

/* namespace reads */
