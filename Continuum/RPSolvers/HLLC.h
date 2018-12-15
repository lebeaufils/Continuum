/*
 HLLC solver for euler equations
 */

#ifndef HLLC_H_
#define HLLC_H_

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

typedef Eigen::Vector3d vector;

class HLLC
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
	HLLC(double, eulerTests); //double c, int N, double L || courant number/ number of cells/ domain length

	void boundary_conditions();
	void initial_conditions(eulerTests, IdealGas);
	void initial_conditions(eulerTests, JWL);
	void solver(IdealGas, eulerTests);
	void solver(JWL, eulerTests);
	void output(IdealGas);
	void output(JWL);
};


#endif /* HLLC_H_ */
