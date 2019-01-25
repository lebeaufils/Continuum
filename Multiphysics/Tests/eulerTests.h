/*
 * test_1.h
 *
 *  Created on: 13 Nov 2018
 *      Author: forte
 */

#ifndef TESTS_EULERTESTS_H_
#define TESTS_EULERTESTS_H_

#include <Eigen/Dense>
typedef Eigen::Vector3d vector;

struct eulerTests
{
	const int N;
	const double L;
	double x0;
	double tstop;

	vector initialL;
	vector initialR;

	eulerTests(double N, double L) : N(N), L(L), x0(0.5), tstop(0), initialL(0, 0, 0), initialR(0, 0, 0) {}

	void test1();
	void test2();
	void test3();
	void test4();
	void test5();
};


#endif /* TESTS_EULERTESTS_H_ */

