#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <Eigen/Dense>
//#include <iostream>
#include "../EOS/EOS.h"

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

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

	StateFunctions* state_function;

	Euler1D() : N(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), state_function(NULL) {}
	//~Euler1D() {}
};






#endif /* VARIABLES_H_ */