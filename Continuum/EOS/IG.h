/*
Ideal Gas EOS
 */

#ifndef EOS_IG_H_
#define EOS_IG_H_

#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;

struct IdealGas
{
	const double y;

	IdealGas() : y(1.4) {}

	double internalE(Eigen::MatrixXd, int);
	double Pressure(Eigen::MatrixXd, int);
	double soundspeed(Eigen::MatrixXd, int);
	double PressureScalar(vector); //used for muscl
	double soundspeedScalar(vector);
	vector f(vector);

};


#endif /* EOS_IG_H_ */
