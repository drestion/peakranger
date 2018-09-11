/*
 * distributions.cpp
 *
 *  Created on: May 6, 2011
 *      Author: xin
 */

#include "distributions.h"

#include <boost/math/distributions/binomial.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>

distributions::distributions() {


}

double distributions::binomCDF(uint32_t n, uint32_t k, double prob) {
	using boost::math::binomial;
	binomial binom(n, prob);
	return cdf(binom, k);
}

double distributions::binomPDF(uint32_t n, uint32_t k, double prob) {
	using boost::math::binomial;
	binomial binom(n, prob);
	return pdf(binom, k);
}
