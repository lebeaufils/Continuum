/*
Stiffened Gas EOS
 */

#ifndef EOS_SG_H_
#define EOS_SG_H_

#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;

struct StiffenedGas
{
	double y;
	double Pref;

	StiffenedGas() : y(1.4) {}

	void GetGamma();
	double internalE(Eigen::MatrixXd, int);
	double Pressure(Eigen::MatrixXd, int);
	double soundspeed(Eigen::MatrixXd, int);
	double PressureScalar(vector); //used for muscl
	double soundspeedScalar(vector);
	vector f(vector);

};


#endif /* EOS_SG_H_ */