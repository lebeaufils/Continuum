#ifndef FDSCHEMES_H_
#define FDSCHEMES_H_

#include <iostream>
#include <vector>
#include <sstream>

class FiniteDifference{
protected:
	//number of cells
	int N;

	//Domain parameters
	double L; //length
	double a; // "Velocity"
	double c; //CFL conditions

	//Cell width
	//double dx;
	//double dt;

	//Arrays
	std::vector<double> X;
	std::vector<double> u;
	std::vector<double> u_1; //time-step n+1

public:
	//Cell width
	double dx;
	double dt;

	virtual void initial_conditions_step();
	virtual void initial_conditions_square();
	virtual void boundary_conditions();
	virtual void output(std::string);
	virtual void plotname(double t);

	FiniteDifference(int, double, double, double);
	//virtual ~FiniteDifference(); //{ std::cout << "Destructor of FD" << std::endl; }
};

class LxF : public virtual FiniteDifference
{
public:
	LxF(int, double, double, double);
	void lax_friedrichs();
	virtual void plotname(double t);
};

class LxW : public virtual FiniteDifference
{
public:
	LxW(int, double, double, double);
	void lax_wendroff();
	virtual void plotname(double t);
};

class WB : public virtual FiniteDifference
{
public:
	WB(int, double, double, double);
	void warming_beam();
	virtual void initial_conditions_step();
	virtual void output(std::string);
	virtual void boundary_conditions();
	virtual void plotname(double t);

	//virtual ~WB(); //std::cout << "Destructor of WB" << std::endl;}
};

#endif /* FDSCHEMES_H_ */
