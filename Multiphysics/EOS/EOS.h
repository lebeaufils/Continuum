/*
Ideal Gas EOS
 */

#ifndef EOS_IG_H_
#define EOS_IG_H_

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;

struct EOS
{
	double y;

	EOS() : y(1.4) {}
	virtual ~EOS() {};

	virtual void GetGamma() = 0;

	double internalE(Eigen::MatrixXd, int);
	double soundspeed(Eigen::MatrixXd, int);
	double soundspeedScalar(vector);
	vector f(vector);
	virtual double Pressure(Eigen::MatrixXd, int) = 0;
	virtual double PressureScalar(vector) = 0; //used for muscl
	virtual vector conservedVar(vector) = 0;
};

struct IdealGas : public virtual EOS
{
	IdealGas() : EOS() {}
	~IdealGas() {};

	void GetGamma();

	double Pressure(Eigen::MatrixXd, int);
	double PressureScalar(vector); //used for muscl
	vector conservedVar(vector);
};

struct StiffenedGas : public virtual EOS
{
	double Pref;

	StiffenedGas() : EOS(), Pref(0) {}
	~StiffenedGas() {};

	void GetGamma();

	double Pressure(Eigen::MatrixXd, int);
	double PressureScalar(vector); //used for muscl
	vector conservedVar(vector);
};


#endif /* EOS_IG_H_ */
