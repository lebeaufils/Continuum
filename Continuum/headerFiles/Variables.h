#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <Eigen/Dense>
#include <vector>
//#include <iostream>
#include "EOS.h"

typedef Eigen::Vector3d vector;
typedef Eigen::Vector4d vector4;
typedef Eigen::MatrixXd matrix;
typedef Eigen::Matrix<vector4, Eigen::Dynamic, Eigen::Dynamic> vecarray;

//-----------------------------------------------------------------------------------------------
//	Data storage
//-----------------------------------------------------------------------------------------------
struct Domain1D
{
	int N;
	double L;
	double dt;
	double dx;
	double tstop;

	matrix X;

	Domain1D() : N(0), L(0), dt(0), dx(0), tstop(0), X(0, 0) {}
	Domain1D(int N) : N(N), L(1.0), dt(0), dx(0), tstop(0), X(0, 0) {}
};

struct Domain2D
{
	int Nx;
	int Ny;
	double Lx;
	double Ly;
	double dt;
	double dx;
	double dy;
	double tstop;

	matrix X;

	Domain2D() : Nx(0), Ny(0), Lx(0), Ly(0), dt(0), dx(0), tstop(0), X(0, 0) {}
	Domain2D(int N) : Nx(N), Ny(N), Lx(1.0), Ly(1.0), dt(0), dx(0), tstop(0), X(N, N) {}
	Domain2D(int Nx, int Ny) : Nx(Nx), Ny(Ny), Lx(1.0), Ly(1.0), dt(0), dx(0), tstop(0), X(Nx, Ny) {}
};

struct Euler1D
{
	//domain parameters
	//int N;
	//double dt;
	//double dx;
	//double x0;
	//double tstop;

	//Variable Matrix
	//matrix X; //Domain
	matrix U; //conserved variables
	matrix F; //flux

	std::shared_ptr<StateFunctions> state_function;

	//Euler1D() : N(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), state_function(NULL) {}
	Euler1D() : U(0, 0), F(0, 0), state_function(NULL) {}
	~Euler1D() {}
};


struct Euler2D
{
	//domain parameters
	//int Nx;
	//int Ny;
	//double dt;
	//double dx;
	//double dy;
	//double x0;
	//double tstop;

	//Variable Matrix
	//matrix X; //Domain
	vecarray U; //conserved variables
	vecarray F; //flux for the x derivative
	vecarray G; // "" y derivative

	std::shared_ptr<StateFunctions> state_function;

	//Euler2D() : Nx(0), Ny(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), G(0, 0), state_function(NULL) {}
	Euler2D() : U(0, 0), F(0, 0), G(0, 0), state_function(NULL) {}
	~Euler2D() {}

	//functions to return row i
	template <typename T>
	T get_row(T, int);
	template <typename T>
	T get_column(T, int);
	template<typename T>
	void display(T);
	template<typename T>
	T swap_xy(T); //swaps the order of velocity in x-y directions
};

struct RB_2D
{
	//collection of level sets, one for each object
	//Eigen::Array<matrix, 1, Eigen::Dynamic> levelset_array;
	std::vector<matrix> levelset_array;
	Euler2D fluid;
	Euler2D rigidbody;

	RB_2D() : levelset_array(0), fluid(), rigidbody() {}
};

struct Polygon
{	
	//for an n-sided polygon
	int n;
	//Eigen::Array<Eigen::Array<double,1,2>, Eigen::Dynamic, 1> vertices;
	//Eigen::Array<Eigen::Array<int,1,2>, Eigen::Dynamic, 1> edges;
	std::vector<Eigen::Array<double,1,2>, Eigen::aligned_allocator<Eigen::Array<double,1,2>>> vertices;
	std::vector<Eigen::Array<int,1,2>, Eigen::aligned_allocator<Eigen::Array<int,1,2>>> edges;
	//faces for 3D

	Polygon(int n) : n(n), vertices(n), edges(n) {}

};


#include "Variables.tcc"




#endif /* VARIABLES_H_ */