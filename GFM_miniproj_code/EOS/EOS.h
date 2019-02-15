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
	Eigen::Matrix<double, 14, 1> C;

	EOS();
	virtual ~EOS() {};

	virtual void GetGamma() = 0;

	double internalE(Eigen::MatrixXd, int);
	double soundspeed(Eigen::MatrixXd, int);
	double soundspeedScalar(vector);
	vector f(vector);
	virtual double Pressure(Eigen::MatrixXd, int) = 0;
	virtual double PressureScalar(vector) = 0; //used for muscl
	virtual vector conservedVar(vector) = 0;

	void testing();

	//Functions for exact solver
	virtual void y_constants(vector) = 0;
	virtual double fk(double) = 0;
	virtual double fprimek(double) = 0;
	//virtual double f(double, vector) = 0;
	//virtual double fprime(double, vector) = 0;
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
	double fk(double);
	double fprimek(double);
	//double f(double, vector, vector);
	//double fprime(double, vector, vector);

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
	double fk(double);
	double fprimek(double);
	//double f(double, vector, vector);
	//double fprime(double, vector, vector);
};


#endif /* EOS_IG_H_ */
