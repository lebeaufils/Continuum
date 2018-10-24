#ifndef FVSCHEMES_H_
#define FVSCHEMES_H_

#include <iostream>
#include <vector>
#include <sstream>

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
	std::vector<double> u;
	std::vector<double> u_1; //time-step n+1
	std::vector<double> f; //flux

public:
	//Cell width
	double dx;
	double dt;

	virtual void initial_conditions_step();
	//virtual void initial_conditions_square();
	virtual void boundary_conditions();
	virtual void output(std::string);
	virtual void plotname(double t);

	FiniteVolume(int, double, double, double);
	//virtual ~FiniteVolume();
};

//centered scheme
class FORCE : public virtual FiniteVolume
{
public:
	FORCE(int, double, double, double);
	void solver();
	//virtual ~FORCE();
};

class linearadvection : public virtual FORCE
{
public:
	linearadvection(int, double, double, double);
	double flux(double);
	void evaluate_flux();
	//virtual void plotname(double t);
	//virtual ~linearadvection();
};


#endif /* FVSCHEMES_H_ */
