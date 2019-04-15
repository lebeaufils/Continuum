/*
 * test_1.h
 *
 *  Created on: 13 Nov 2018
 *      Author: forte
 */

#ifndef TESTS_EULERTESTS_H_
#define TESTS_EULERTESTS_H_

#include <Eigen/Dense>
#include "EOS.h"
#include "Variables.h"

typedef Eigen::Vector3d vector;
typedef Eigen::Vector4d vector4;
typedef Eigen::MatrixXd matrix;



//-----------------------------------------------------------------------------------------------
//	Tests
//-----------------------------------------------------------------------------------------------

struct standardTests
{
	//int N;
	//double L;
	//double x0;
	//double tstop;

	matrix initialL;
	matrix initialR;

	//standardTests(int N) : N(N), L(1.0), x0(0.5), tstop(0), initialL(3, 1), initialR(3, 1) {}
	standardTests() : initialL(3, 1), initialR(3, 1) {}
	virtual ~standardTests() {}

	int Get_Switch();
	void switch_resolution(Domain1D&);
};

struct eulerTests : public virtual standardTests
{
	Domain1D domain;
	Euler1D var;

	double x0; //interface location
	//only a single interface is supported for standard toro tests

	eulerTests(int N) : standardTests(), domain(N), var() {}
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

struct eulerTests2D : public virtual standardTests
{
	Domain2D domain;
	Euler2D var;

	//Additional domain parameters
		//If x and y have different lengths/ number of cells
		//int Ny;
		//double Ly;
		//Interface location as a function of x and y (a list of interfacial cells)
		Eigen::Array<bool,Eigen::Dynamic,Eigen::Dynamic> interface;

	//eulerTests2D(int N) : standardTests(N), var(), Ny(N), Ly(L), interface(N, Ny) 
	eulerTests2D(int N) : standardTests(), domain(N), var(), interface(N, N) {
		initialL.resize(4, 1);
		initialR.resize(4, 1);
	}
	eulerTests2D(int Nx, int Ny) : standardTests(), domain(Nx, Ny), var(), interface(Nx, Ny) {
		initialL.resize(4, 1);
		initialR.resize(4, 1);
	}
	virtual ~eulerTests2D() {}

	void test1();
	void test2();
	void test3();
	void test4(); //cylindrical shock test

	//void switch_test();	
};

struct rigidTests : public virtual standardTests
{
	Domain1D domain;
	Euler1D fluid;
	Euler1D rigidbody;

	matrix interfacelist; //list of interface coordinates

	rigidTests(int N) : standardTests(), domain(N), fluid(), rigidbody(), interfacelist(0, 0) {}
	virtual ~rigidTests() {}

	void test1(); 
};


/*struct gfmTests : public virtual standardTests
{
	//For now, only capable of dealing with 2 different materials.
	Euler1D var1;
	Euler1D var2;

	int number_of_discontinuities;
	bool stiffgas = false; //flag for stiffened gas

	//double yL;
	//double yR;
	//double yM1;
	//double yM2;
	//double Pref1;
	//double Pref2; //these values are stored inside the variable struce

	double x1;
	double x2;

	//storage of initial conditions for up to 4 discontinuities.
	vector initialM1;
	vector initialM2;

	gfmTests() : standardTests(), number_of_materials(2),
	x1(0), x2(0), initialM1(0, 0, 0), initialM2(0, 0, 0)  {}
	gfmTests(double N, double L) : standardTests(N, L), number_of_materials(2),
	x1(0), x2(0), initialM1(0, 0, 0), initialM2(0, 0, 0) {}

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

