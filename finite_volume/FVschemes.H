#ifndef FVSCHEMES_H_
#define FVSCHEMES_H_

#include <iostream>
#include <vector>
#include <sstream>
//#include <Eigen/Dense> //dense matrices only. use sparse if needed

class FiniteVolume{
protected:
	//number of cells
	int N;

	//Domain parameters
	double L; //length
	double a; // "Velocity"
	double c; //CFL conditions

	//Arrays
	std::vector<double> X;

public:
	//Cell width
	double dx;
	double dt;

	FiniteVolume(int, double, double, double);
	//virtual ~FiniteVolume();
};

class scalar : public virtual FiniteVolume
{
protected:
	std::vector<double> u;
	std::vector<double> u_1; //time-step n+1
	std::vector<double> f; //flux

public:
	scalar(int, double, double, double);

	virtual void initial_conditions_step();
	virtual void initial_conditions_square();
	virtual void boundary_conditions();
	virtual void output(std::string);
	virtual void plotname(double t);
};

class euler : public virtual FiniteVolume
{
	/*Eigen::MatrixXd U;
	Eigen::MatrixXd U_1;
	Eigen::MatrixXd F;*/
	std::vector<std::vector<double> > U;
	std::vector<std::vector<double> > U_1;
	std::vector<std::vector<double> > F;

public:
	euler(int, double, double, double);

	virtual void initial_conditions_step();
	virtual void boundary_conditions();
	virtual void output(std::string);
	virtual void plotname(double t);
};

//centered scheme
class scalarFORCE : public virtual scalar
{
public:
	scalarFORCE(int, double, double, double);
	void solver();
	//virtual ~FORCE();
};

class eulerFORCE: public virtual euler
{
public:
	eulerFORCE(int, double, double, double);

};

class linearadvection : public virtual scalarFORCE
{
public:
	linearadvection(int, double, double, double);
	double flux(double);
	void evaluate_flux();
	virtual void plotname(double t);
	//virtual ~linearadvection();
};

/*class Euler: public virtual eulerFORCE
{
	std::vector<double> p; //density
	std::vector<double> pu; //density*velocity
	std::vector<double> E; //energy
};*/


#endif /* FVSCHEMES_H_ */
