/*
 * test_1.h
 *
 *  Created on: 13 Nov 2018
 *      Author: forte
 */

#ifndef TESTS_EULERTESTS_H_
#define TESTS_EULERTESTS_H_

#include <Eigen/Dense>
#include "../EOS/EOS.h"
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
	void test1_stationary();
	void test2();
	void test3();
	void test4();
	void test5();
	void test6();
	void test7();

	//custom tests
	//void testN(); //user input
	//void testS(); //settings file
};


struct gfmTests : public virtual standardTests
{
	int number_of_materials;
	double yL;
	double yR;
	double yM1;
	double yM2;
	double Pref1;
	double Pref2;

	double x1;
	double x2;
	vector initialM1;
	vector initialM2;

	gfmTests() : standardTests(0, 0), number_of_materials(2), yL(1.4), yR(1.4), yM1(0), yM2(0), Pref1(0), Pref2(0),
	x1(0), x2(0), initialM1(0, 0, 0), initialM2(0, 0, 0)  {}
	gfmTests(double N, double L) : standardTests(N, L), number_of_materials(2), yL(1.4), yR(1.4), yM1(0), yM2(0), Pref1(0), Pref2(0),
	x1(0), x2(0), initialM1(0, 0, 0), initialM2(0, 0, 0) {}

	void set_number_of_cells(int); //changing the resolution. 

	void test_example_1();
	void testA();
	void testB();
	void testC();
	void testB_Wang();
	void testD();
	void testE();
	void testF();

	//custom tests
	//void testN(); //user input
	//void testS(); //settings file

	//Set EOS parameters
	void set_EOS(EOS*, EOS*);
};
	



#endif /* TESTS_EULERTESTS_H_ */

