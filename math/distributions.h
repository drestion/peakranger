/*
 * distributions.h
 *
 *  Created on: May 6, 2011
 *      Author: xin
 */

#ifndef DISTRIBUTIONS_H_
#define DISTRIBUTIONS_H_
#include <stdint.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime>
#include <stdlib.h>
#include <vector>

namespace {
boost::mt19937 gen(std::time(0));
}

class distributions {
public:
    distributions();

    static double binomCDF(uint32_t n, uint32_t k, double prob);
    static double binomPDF(uint32_t n, uint32_t k, double prob);

    /*
     * It seems that in real world it is faster
     * than the std version
     */
    static void random_int_boost(uint32_t l, uint32_t h, uint32_t size,
            std::vector<uint32_t>& result) {
        //boost::mt19937 gen(std::time(0));
        boost::mt19937 gen(12345);
        for (uint32_t i = 0; i < size; ++i) {
            boost::uniform_int<> dist(l, h);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rnd(
                    gen, dist);
            uint32_t r = (uint32_t) (rnd());
            result.push_back(r);
        }

    }

    /*
     * The divide & conquer version of boost,slower than boost
     */
    static void random_int_boost_dc(uint32_t l, uint32_t h, uint32_t size,
            std::vector<uint32_t>& result) {

        if (size == 1) {

            boost::uniform_int<> dist(l, h);
            boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rnd(
                    gen, dist);
            uint32_t r = (uint32_t) (rnd());
            result.push_back(r);
            return;
        }

        uint32_t half = size / 2;

        // Run the function recursively on both halves.
        random_int_boost_dc(l, h, half, result);

        random_int_boost_dc(l, h, size - half, result);

    }

    static void random_int_std(uint32_t l, uint32_t h, uint32_t size,
            std::vector<uint32_t>& result) {
//        srand((int) time(NULL));
    	srand(12345);
        for (uint32_t i = 0; i < size; ++i) {
            uint32_t k;
            double d;
            d = (double) rand() / ((double) RAND_MAX + 1);
            k = (uint32_t) (d * (h - l + 1));
            uint32_t r = l + k;
            result.push_back(r);
        }
    }

};

#endif /* DISTRIBUTIONS_H_ */
