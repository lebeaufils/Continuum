
#ifndef SOLVERS_H_
#define SOLVERS_H_

//Equation of states and tests
//#include "../EOS/JWL.h"
#include "../EOS/EOS.h"
#include "../Tests/eulerTests.h"
//#include "Variables.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <map>

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

//--------------------------------------------------------------------------------
//	MUSCL
//--------------------------------------------------------------------------------
enum slopeLimiter {MinBee, VanLeer, SuperBee, Quit};

struct MUSCL
{
	static void boundary_conditions(Euler1D&);
	static void initial_conditions_1D(eulerTests&);

	//-----Slope limiters-----
	static vector superBee(matrix, int);
	static vector vanLeer(matrix, int);
	static vector minBee(matrix, int);
	//------------------------
	//-----Data Reconstruction-----
	static slopeLimiter getLimiter();
	static void data_reconstruction(matrix, slopeLimiter, matrix&, matrix&, int);
	//-----------------------------

	static void compute_fluxes(Euler1D&, int, matrix, matrix, double&);
	static void conservative_update_formula(Euler1D&, int);

	static void solver(Euler1D&, double);
	static void output(Euler1D&);
};

struct EXACT 
{
//------------------------------------------------
//	Input Parameters
//------------------------------------------------
	int N;
	double dx;
	double x0;

	//VDomain and initial conditions
	//matrix X; //Domain
	matrix W;
	vector WL;
	vector WR;

//------------------------------------------------
//	Exact solver 
//------------------------------------------------
	//Tolerance level (Mixed error testing)
	double TOL;
	//constants
	double y;
	double cL, cR; //soundspeed of left and right states
	double CONST1; // y-1 / 2y
	double CONST2; // y+1 / 2y
	double CONST3; // 2*y / y-1
	double CONST4; // 2 / y-1
	double CONST5; // 2 / y+1
	double CONST6; // y-1 / y+1
	double CONST7; // y-1 / 2
	double CONST8; // y-1

public:
	EXACT(eulerTests Test) : N(1000), dx(Test.L/N), x0(Test.x0), W(N, 3), TOL(1e-6), 
		y(Test.var.state_function->y), cL(0), cR(0), CONST1(0),CONST2(0), CONST3(0), CONST4(0), CONST5(0), CONST6(0), CONST7(0), CONST8(0) {
			WL = Test.initialL;
			WR = Test.initialR;
		}
	
	//EXACT(int, double, double, double, vector, vector); //unfinished
	//EXACT(double, eulerTests);
	//virtual ~EXACT() {};

	//compute constants
	void initial_conditions();

	//Riemann Problem: Equations for Pressure and Particle Velocity

	//Pressure positivity condition Toro pg 127
		//Direct evaluation of f(p) gives the pressure positivity condition
	void check_pressure_pos_condition();

	//void fL(double);
	double fk(double, vector);
	double f(double);
	double fkprime(double, vector);
	double fprime(double);
	double newton_raphson(double);
	double relative_pressure_change(double, double);

	double compute_star_pressure(); //Newton-Raphson Method
	double compute_star_velocity(double);
	double compute_shock_density(vector, double);
	double compute_rarefraction_density(vector, double);
	void sampling(double);
	void output();
	void solver(eulerTests&);

	//double compute_star_velocity(vector, vector, EOS*, double);
	//double compute_star_shockdensity_k(vector, EOS*, double);
	//double compute_mass_flux_k(vector, EOS*, double);
	//double compute_left_shock_speed(vector, EOS*, double);
	//double compute_right_shock_speed(vector, EOS*, double);
	//double compute_rarefraction_speed();

	//if P* > PL --> shock wave
	//compute shock speed,
	//W = W*shock if x/t > SL and x/t< Ustar
	//if P* < PL --> rarefractionwave
	//compute speed head and speed tail
	//W = WL, WLfan, WLfan*....

	//c^2 = y*(P+P0)/d


};

#endif /* SOLVERS_H_ */
