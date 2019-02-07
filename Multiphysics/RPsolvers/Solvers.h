/*
 * Solvers.h
 *
 *  Created on: 24 Jan 2019
 *      Author: forte
 */

#ifndef SOLVERS_H_
#define SOLVERS_H_

//Equation of states and tests
//#include "../EOS/JWL.h"
#include "../EOS/EOS.h"
#include "../Tests/eulerTests.h"

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

class RPsolvers
{
protected:
	const double CFL; //Courant number

	//domain parameters
	int N;
	int count;
	double dt;
	double dx;

	//Variable Matrix
	matrix X; //Domain
public:
	matrix U; //Conserved quantities d, du, E
	matrix F; //Flux
	double Smax; //Maximum soundspeed

	RPsolvers() : CFL(0), N(0), count(0), dt(0), dx(0), X(0, 0), U(0, 0), F(0, 0), Smax(0){}
	RPsolvers(gfmTests,  int, int); //without CFL, for gfm
	RPsolvers(double, eulerTests,  int, int); //double c, int N, double L || courant number/ number of cells/ domain length
	virtual ~RPsolvers() {};

	void conservative_update_formula(int);
	void conservative_update_formula(double, double, int);
	//Abstract classes
	virtual void boundary_conditions() = 0;
	virtual void initial_conditions(EOS*, eulerTests) = 0;
	//virtual void initial_conditions(JWL, eulerTests) = 0;
	virtual void compute_fluxes(EOS*, int) = 0;
	//virtual void compute_fluxes(JWL, eulerTests) = 0;
	virtual void solver(EOS*, eulerTests) = 0;
	virtual void output(EOS*) = 0;
	//virtual void output(JWL) = 0;
};

/*
class FiniteVolume : public virtual RPsolvers
{	
protected:
	matrix X; //Domain
public:
	matrix U; //Conserved quantities d, du, E
	matrix F; //Flux

	FiniteVolume() : RPsolvers(){}
	FiniteVolume(gfmTests,  int, int);
	FiniteVolume(double, eulerTests,  int, int);
	virtual ~FiniteVolume() {};

	virtual void conservative_update_formula(int);
	virtual void conservative_update_formula(double, double, int);
	//Abstract classes
	virtual void boundary_conditions() = 0;
	virtual void initial_conditions(EOS*, eulerTests) = 0;
	virtual void compute_fluxes(EOS*, int) = 0;
	virtual void solver(EOS*, eulerTests) = 0;
	virtual void output(EOS*) = 0;
};
*/

/*--------------------------------------------------------------------------------
 * HLLC
 --------------------------------------------------------------------------------*/

class HLLC : public virtual RPsolvers
{
public:
	HLLC(gfmTests);
	HLLC(double, eulerTests);
	virtual ~HLLC() {};

	virtual void boundary_conditions();
	virtual void initial_conditions(EOS*, eulerTests);
	//virtual void initial_conditions(JWL, eulerTests);
	virtual void compute_fluxes(EOS*, int);
	//virtual void compute_fluxes(JWL, eulerTests);
	virtual void solver(EOS*, eulerTests);
	virtual void output(EOS*);
	//virtual void output(JWL);
};

/*--------------------------------------------------------------------------------
 * MUSCL
 --------------------------------------------------------------------------------*/
enum slopeLimiter {MinBee, VanLeer, SuperBee, Quit};

class MUSCL : public virtual RPsolvers
{
	matrix ULi; 
	matrix URi;

public:
	MUSCL(gfmTests);
	MUSCL(double, eulerTests);
	virtual ~MUSCL() {};

	virtual void boundary_conditions();
	virtual void initial_conditions(EOS*, eulerTests);
	//virtual void initial_conditions(JWL, eulerTests);

	//-----Slope limiters-----
	vector superBee(int);
	vector vanLeer(int);
	vector minBee(int);
	//------------------------
	//-----Data Reconstruction-----
	slopeLimiter getLimiter();
	void data_reconstruction(slopeLimiter);
	//-----------------------------

	/*
	vector f(vector, EOS*);
	vector f(vector, JWL);
	*/
	virtual void compute_fluxes(EOS*, int);
	//virtual void compute_fluxes(JWL, eulerTests);
	virtual void solver(EOS*, eulerTests);
	virtual void output(EOS*);
	//virtual void output(JWL);

};

class EXACT //: public virtual RPsolvers
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
	EXACT(eulerTests Test, EOS* IG) : N(Test.N), dx(Test.L/Test.N), x0(Test.x0), W(N, 3), TOL(1e-6), 
		y(IG->y), cL(0), cR(0), CONST1(0),CONST2(0), CONST3(0), CONST4(0), CONST5(0), CONST6(0), CONST7(0), CONST8(0) {
			WL = Test.initialL;
			WR = Test.initialR;
		}
	
	EXACT(int, double, double, double, vector, vector); //unfinished
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
