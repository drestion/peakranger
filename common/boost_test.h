/*
 * boost_test.h
 *
 *  Created on: May 4, 2012
 *      Author: xin
 */

#ifndef BOOST_TEST_H_
#define BOOST_TEST_H_


#define Test(t) BOOST_AUTO_TEST_CASE(t)
#define SuiteBegin(t,F) BOOST_FIXTURE_TEST_SUITE(T,F)
#define SuiteEnd BOOST_AUTO_TEST_SUITE_END
#define checkEqual BOOST_CHECK_EQUAL
#define Fail BOOST_FAIL


#endif /* BOOST_TEST_H_ */
