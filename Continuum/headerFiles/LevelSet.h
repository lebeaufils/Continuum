#ifndef LEVELSET_H_
#define LEVELSET_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

#include "eulerTests.h"
#include "Variables.h"

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

class LevelSetFunction1D
{
protected:
	//domain parameters
	/*int N;
	matrix X;
	double dx;
	double x0;
	double x1;
	double x2;*/

	matrix phi;
	//int sgn; //sign of the levelset

public:
	//LevelSetFunction(gfmTests);
	LevelSetFunction();

	int get_sgn(double);
	void boundary_conditions(Euler1D);
	void signed_distance_function_1D(int);
	void signed_distance_function_1D_2(int); //2 discontinuties
	void signed_distance_function_1D_3(int); //3 discontinuities
	
	//first order upwind hamilton-jacobi method
	double HJ_FirstOrder(double, double, int); //velocity
	//void HJ_ENO(VectorN2, double);
	//void HJ_WENO(VectorN2, double);

	//reinitialisation
	void reinitialisation();
	//void reconstruction();
};

class LevelSetFunction2D
{
protected:
	//domain parameterz
	matrix phi; //List of level set values

public:
	LevelSetFunction2D();

	int get_sgn(double);
	void boundary_conditions();
	void signed_distance_function_2D(int);

	//first order upwind hamilton-jacobi method
	double HJ_FirstOrder(double, double, int); //velocity
	//void HJ_ENO(VectorN2, double);
	//void HJ_WENO(VectorN2, double);

	//reinitialisation
	void reinitialisation(); //fast marching method
	//void reconstruction();
};


#endif /* LEVELSET_H_ */