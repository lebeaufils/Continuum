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
	int N;
	int count;
	double dt;
	double dx;

	//Variable Matrix
	matrix X; //Domain
	matrix W;

	double TOL;

public:
	EXACT(eulerTests Test) : N(Test.N), count(0), dt(0), dx(Test.L/Test.N), X(N, 1), W(N, 3), TOL(1e-6) {}
	//EXACT(double, eulerTests);
	//virtual ~EXACT() {};

	void initial_conditions(eulerTests);

	//Riemann Problem: Equations for Pressure and Particle Velocity
	//void fL(double);
	double fk(double, vector, EOS*);
	double f(double, vector, vector, EOS*);
	double fkprime(double, vector, EOS*);
	double fprime(double, vector, vector, EOS*);
	double newton_raphson(double, vector, vector, EOS*);
	double compute_star_pressure(vector, vector, EOS*); //Newton-Raphson Method
	double relative_pressure_change(double, double);




};

#endif /* SOLVERS_H_ */
