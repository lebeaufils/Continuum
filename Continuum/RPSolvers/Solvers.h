/*
 * Solvers.h
 *
 *  Created on: 24 Jan 2019
 *      Author: forte
 */

#ifndef SOLVERS_H_
#define SOLVERS_H_

//Equation of states and tests
#include "../EOS/JWL.h"
#include "../EOS/IG.h"
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
	Eigen::MatrixXd X; //Domain
	Eigen::MatrixXd U; //Conserved quantities d, du, E
	Eigen::MatrixXd F; //Flux

public:
	RPsolvers(double, eulerTests, int, int); //double c, int N, double L || courant number/ number of cells/ domain length
	virtual ~RPsolvers() {};
	//Abstract classes
	virtual void boundary_conditions() = 0;
	virtual void initial_conditions(IdealGas, eulerTests) = 0;
	virtual void initial_conditions(JWL, eulerTests) = 0;
	virtual void solver(IdealGas, eulerTests) = 0;
	virtual void solver(JWL, eulerTests) = 0;
	virtual void output(IdealGas) = 0;
	virtual void output(JWL) = 0;
};

/*--------------------------------------------------------------------------------
 * HLLC
 --------------------------------------------------------------------------------*/

class HLLC : public virtual RPsolvers
{
public:
	HLLC(double, eulerTests);
	~HLLC() {};

	virtual void boundary_conditions();
	virtual void initial_conditions(IdealGas, eulerTests);
	virtual void initial_conditions(JWL, eulerTests);
	virtual void solver(IdealGas, eulerTests);
	virtual void solver(JWL, eulerTests);
	virtual void output(IdealGas);
	virtual void output(JWL);
};

/*--------------------------------------------------------------------------------
 * MUSCL
 --------------------------------------------------------------------------------*/

class MUSCL : public virtual RPsolvers
{
public:
	MUSCL(double, eulerTests);
	~MUSCL() {};

	virtual void boundary_conditions();
	virtual void initial_conditions(IdealGas, eulerTests);
	virtual void initial_conditions(JWL, eulerTests);

	//-----Slope limiters-----
	vector superBee(int);
	vector vanLeer(int);
	vector minBee(int);
	//------------------------
	/*
	vector f(vector, IdealGas);
	vector f(vector, JWL);
	*/
	virtual void solver(IdealGas, eulerTests);
	virtual void solver(JWL, eulerTests);
	virtual void output(IdealGas);
	virtual void output(JWL);

};


#endif /* SOLVERS_H_ */
