#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <Eigen/Dense>
//#include <iostream>
#include "EOS.h"

typedef Eigen::Vector3d vector;
typedef Eigen::Vector4d vector4;
typedef Eigen::MatrixXd matrix;
typedef Eigen::Matrix<vector4, Eigen::Dynamic, Eigen::Dynamic> vecarray;

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


struct Euler2D
{
	//domain parameters
	int Nx;
	int Ny;
	double dt;
	double dx;
	double dy;
	double x0;
	double tstop;

	//Variable Matrix
	matrix X; //Domain
	vecarray U; //conserved variables
	vecarray F; //flux for the x derivative
	vecarray G; //

	std::shared_ptr<StateFunctions> state_function;

	Euler2D() : Nx(0), Ny(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), G(0, 0), state_function(NULL) {}
	~Euler2D() {}

	//functions to return row i
	template <typename T>
	T get_row(T, int);
	template <typename T>
	T get_column(T, int);
	template<typename T>
	void display(T);
};

#include "Variables.tcc"




#endif /* VARIABLES_H_ */