/*
 * MUSCL.h
 *
 *  Created on: 15 Dec 2018
 *      Author: forte
 */

#ifndef MUSCL_H_
#define MUSCL_H_

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

class MUSCL
{
	const double CFL; //Courant number
	Eigen::MatrixXd X; //Domain
	Eigen::MatrixXd U; //Conserved quantities d, du, E
	Eigen::MatrixXd F; //Flux

	//domain parameters
	int count;
	int N;
	double dx;
	double dt;

public:
	MUSCL(double, eulerTests); //double c, int N, double L || courant number/ number of cells/ domain length

	void boundary_conditions();
	void initial_conditions(IdealGas, eulerTests);
	void initial_conditions(JWL, eulerTests);

	//-----Slope limiters-----
	vector superBee(int);
	vector vanLeer(int);
	vector minBee(int);
	//------------------------

	vector f(vector, IdealGas);
	vector f(vector, JWL);
	void solver(IdealGas, eulerTests);
	void solver(JWL, eulerTests);
	void output(IdealGas);
	void output(JWL);
};





#endif /* MUSCL_H_ */
