#ifndef LEVELSET_H_
#define LEVELSET_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

class LevelSetFunction
{
	//domain parameters
	int N;
	double L;
	double dx;

	//LevelSetFunction
	vector phi;
	vector X;
	int sgn; //sign of the levelset

public:
	LevelSetFunction(eulerTests);

	void boudnary_conditions();
	void initial_conditions();
	void get_sgn();
	/*
	//sign function
	int sgn(double x){
	return (x > 0) - (x < 0);
	Here's a more readable way to do it:

	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;}
	*/
	

	//first order upwind hamilton-jacobi method
	void HJ_FirstOrder(double); //velocity
	void HJ_FirstOrder(double, double);	
	//void HJ_ENO(VectorN2, double);
	//void HJ_WENO(VectorN2, double);

	//reinitialisation
	//void FastMarchingMethod();
};


#endif /* LEVELSET_H_ */