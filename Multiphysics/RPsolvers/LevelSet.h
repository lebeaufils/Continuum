#ifndef LEVELSET_H_
#define LEVELSET_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

#include "../Tests/eulerTests.h"

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

class LevelSetFunction
{
protected:
	//domain parameters
	int N;
	matrix X;
	double dx;
	double x0;
	double x1;
	double x2;

	matrix phi;
	int sgn; //sign of the levelset

public:
	LevelSetFunction(gfmTests);

	int get_sgn(double);
	void boundary_conditions();
	void signed_distance_function_1D(int);
	void signed_distance_function_1D_2(int); //2 discontinuties
	
	//first order upwind hamilton-jacobi method
	double HJ_FirstOrder(double, double, int); //velocity
	//void HJ_ENO(VectorN2, double);
	//void HJ_WENO(VectorN2, double);

	//reinitialisation
	//void FastMarchingMethod();
};


#endif /* LEVELSET_H_ */