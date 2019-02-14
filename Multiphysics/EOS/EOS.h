/*
Ideal Gas EOS
 */

#ifndef EOS_IG_H_
#define EOS_IG_H_

#include <iostream>
#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

struct EOS
{
	double y;
	matrix C;

	EOS() : y(1.4), C(0, 0) {}
	virtual ~EOS() {};

	virtual void GetGamma() = 0;

	double internalE(Eigen::MatrixXd, int);
	double soundspeed(Eigen::MatrixXd, int);
	double soundspeedScalar(vector);
	vector f(vector);
	virtual double Pressure(Eigen::MatrixXd, int) = 0;
	virtual double PressureScalar(vector) = 0; //used for muscl
	virtual vector conservedVar(vector) = 0;
	//for exact soln
	virtual void y_constants(vector) = 0;
};

struct IdealGas : public virtual EOS
{
	IdealGas() : EOS() {}
	~IdealGas() {};

	void GetGamma();

	double Pressure(Eigen::MatrixXd, int);
	double PressureScalar(vector); //used for muscl
	vector conservedVar(vector);

	//exact
	void y_constants(vector);
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

	//exact
	void y_constants(vector);
};


#endif /* EOS_IG_H_ */
