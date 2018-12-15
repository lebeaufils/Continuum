/*
JWL EOS in Mie Gruneisen form
 */

#ifndef EOS_JWL_H_
#define EOS_JWL_H_

#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;

struct JWL
{
	const double Gruneisen;
	const double d0;
	const double A;
	const double B;
	const double R1;
	const double R2;

	JWL() :
		Gruneisen(0.25),
		d0(1840),
		A(854.5e9),
		B(20.5e9),
		R1(4.6),
		R2(1.35)
		{}

	JWL(double g, double d, double a, double b, double r1, double r2) :
		Gruneisen(g),
		d0(d),
		A(a),
		B(b),
		R1(r1),
		R2(r2)
		{}

	double Pref(double);
	double eref(double);
	double internalE(Eigen::MatrixXd, int);
	double Pressure(Eigen::MatrixXd, int);
	double soundspeed(Eigen::MatrixXd, int);
	double PressureScalar(vector);
	double soundspeedScalar(vector);
};


#endif /* EOS_JWL_H_ */
