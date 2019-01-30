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

struct standardTests
{
	int N;
	double L;
	double x0;
	double tstop;

	vector initialL;
	vector initialR;

	standardTests(double N, double L) : N(N), L(L), x0(0.5), tstop(0), initialL(0, 0, 0), initialR(0, 0, 0) {}
};

struct eulerTests : public virtual standardTests
{

	eulerTests(double N, double L) : standardTests(N, L) {}

	void test1();
	void test2();
	void test3();
	void test4();
	void test5();

	//custom tests
	//void testN(); //user input
	//void testS(); //settings file
};


struct gfmTests : public virtual standardTests
{
	int number_of_materials;
	double y1;
	double y2;
	double y3;
	double y4;

	gfmTests() : standardTests(0, 0), number_of_materials(2), y1(1.4), y2(1.4), y3(0), y4(0) {}
	gfmTests(double N, double L) : standardTests(N, L), number_of_materials(2), y1(1.4), y2(1.4), y3(0), y4(0) {}

	void test_example_1();
	void testA();
	void testA_hires();
	void testB();
	void testC();
	void testD();

	//custom tests
	//void testN(); //user input
	//void testS(); //settings file
};
	



#endif /* TESTS_EULERTESTS_H_ */

