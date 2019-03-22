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
//#include "../RPsolvers/Variables.h"

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

//-----------------------------------------------------------------------------------------------
//	Data storage
//-----------------------------------------------------------------------------------------------

struct Euler1D
{
	//domain parameters
	int N;
	double dt;
	double dx;
	double x0;
	double tstop;

	//Variable Matrix
	matrix X; //Domain
	matrix U; //conserved variables
	matrix F; //flux

	std::shared_ptr<StateFunctions> state_function;

	Euler1D() : N(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), state_function(NULL) {}
	~Euler1D() {}

};

//-----------------------------------------------------------------------------------------------
//	Tests
//-----------------------------------------------------------------------------------------------

struct standardTests
{
	int N;
	double L;
	double x0;
	double tstop;

	vector initialL;
	vector initialR;

	standardTests(double N) : N(N), L(1.0), x0(0.5), tstop(0), initialL(0, 0, 0), initialR(0, 0, 0) {}
	virtual ~standardTests() {}

	int Get_Switch();
	void switch_resolution();
};

struct eulerTests : public virtual standardTests
{
	Euler1D var;

	eulerTests(double N) : standardTests(N), var() {}
	virtual ~eulerTests() {}

	void test1();
	void test1_stationary();
	void test2();
	void test3();
	void test4();
	void test5();

	//custom tests
	//void testN(); //user input
	//void testS(); //settings file

	void switch_test();
};

/*
struct gfmTests : public virtual standardTests
{
	int number_of_materials;
	bool stiffgas = false; //flag for stiffened gas

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
	//single material tests
	void test1();
	void test2();
	void test3();
	void test4();
	void test5();
	void test_example_1();

	//Multimaterial tests
	void testA();
	void testB();
	void testC();
	void testB_Wang();
	void testD();
	void testSG();
	void testMach10();
	void testMach10_2();

	void testcase2();

	//user input test choices
	void switch_resolution();
	void switch_test();

	//custom tests
	//void testN(); //user input
	//void testS(); //settings file

	//Set EOS parameters
	//void set_EOS(EOS*, EOS*);
};
	
*/


#endif /* TESTS_EULERTESTS_H_ */

