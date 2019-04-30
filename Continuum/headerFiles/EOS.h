#ifndef EOS_IG_H_
#define EOS_IG_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include <iostream>
#include <memory>
#include <cmath>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;
typedef Eigen::Vector4d vector4;
typedef Eigen::MatrixXd matrix;
typedef Eigen::Matrix<double, 14, 1> Cmatrix;
//typedef std::shared_ptr<StateFunctions> StateFunctionPtr;


enum EOStype {EOS_IG, EOS_SG};
struct StateFunctions
{
	double y;
	StateFunctions() : y(0) {}
	StateFunctions(double y) : y(y) {} 
	virtual ~StateFunctions() {};

	static std::shared_ptr<StateFunctions> create(EOStype type);

	virtual double Pressure(Eigen::MatrixXd, int) = 0;
	virtual double Pressure(vector) = 0; //used for muscl
	virtual double Pressure(vector4) = 0;
	virtual double soundspeed(Eigen::MatrixXd, int);
	virtual double soundspeed(vector);
	virtual double soundspeed(vector4);
	virtual vector conservedVar(vector) = 0;
	virtual vector4 conservedVar2Dx(vector4) = 0;
	virtual vector4 conservedVar2Dy(vector4) = 0;
	virtual vector4 primitiveVar(vector4) = 0;
	double internalE(Eigen::MatrixXd, int);
	double internalE(vector4);
	vector fluxes(vector); //flux
	//2D
	vector4 fluxes(vector4); //flux in the x or y direction

	//Functions for exact solver
	virtual void y_constants(vector) = 0;
	virtual double fk(double) = 0;
	virtual double fprimek(double) = 0;
	//virtual double f(double, vector) = 0;
	//virtual double fprime(double, vector) = 0;
};

struct IdealGas : public virtual StateFunctions
{
	Cmatrix C; //Constants that are related to the isentropic expansion factor
	IdealGas();
	IdealGas(double y);
	~IdealGas() {};

	double Pressure(matrix, int);
	double Pressure(vector); //used for muscl
	double Pressure(vector4);
	//double soundspeed(matrix, int);
	//double soundspeed(vector);
	vector conservedVar(vector);
	//2D
	vector4 conservedVar2Dx(vector4); 
	vector4 conservedVar2Dy(vector4);
	vector4 primitiveVar(vector4);
	//the conserved variables have been reordered to read
		//Density, momentum in the sweep-direction, energy, momentum in the other direction

	//exact
	void y_constants(vector);
	double fk(double);
	double fprimek(double);
	//double f(double, vector, vector);
	//double fprime(double, vector, vector);
};

struct StiffenedGas : public virtual StateFunctions
{
	double Pref;
	Cmatrix C;
	StiffenedGas();
	StiffenedGas(double y, double Pref);
	~StiffenedGas() {};

	double Pressure(matrix, int);
	double Pressure(vector); //used for muscl
	double Pressure(vector4);
	//double soundspeed(matrix, int);
	//double soundspeed(vector);
	vector conservedVar(vector);
	//2D
	vector4 conservedVar2Dx(vector4);
	vector4 conservedVar2Dy(vector4);
	vector4 primitiveVar(vector4);
	//exact
	void y_constants(vector);
	double fk(double);
	double fprimek(double);
	//double f(double, vector, vector);
	//double fprime(double, vector, vector);
};


/*struct StateFunctions
{
	double y;
	Eigen::Matrix<double, 14, 1> C;

	StateFunctions();
	virtual ~StateFunctions() {};

	virtual void GetGamma() = 0;

	double internalE(Eigen::MatrixXd, int);
	vector f(vector);
	virtual double Pressure(Eigen::MatrixXd, int) = 0;
	virtual double Pressure(vector) = 0; //used for muscl
	virtual double soundspeed(Eigen::MatrixXd, int);
	virtual double soundspeed(vector);
	virtual vector conservedVar(vector) = 0;

	void testing();

	//Functions for exact solver
	virtual void y_constants(vector) = 0;
	virtual double fk(double) = 0;
	virtual double fprimek(double) = 0;
	//virtual double f(double, vector) = 0;
	//virtual double fprime(double, vector) = 0;
};

struct StateFunctions //computes ideal gas state functions
{
	static double pressure(matrix, int, IdealGas);
	static double pressure(vector, IdealGas); //used for muscl
	static double soundspeed(matrix, int, IdealGas);
	static double soundspeed(vector, IdealGas);
	static vector conserved(vector, IdealGas);

	//exact
	static Cmatrix y_constants(vector, IdealGas);
	static double fk(double, IdealGas);
	static double fprimek(double, IdealGas);

	//double f(double, vector, vector);
	//double fprime(double, vector, vector);

};
*/
/*
struct SGStateFunctions
{
	double Pref;

	StiffenedGas() : StateFunctions(), Pref(0) {}
	~StiffenedGas() {};

	double Pressure(matrix, int);
	double Pressure(vector); //used for muscl
	double soundspeed(matrix, int);
	double soundspeed(vector);
	vector conservedVar(vector);

	//exact
	void y_constants(vector);
	double fk(double);
	double fprimek(double);
	//double f(double, vector, vector);
	//double fprime(double, vector, vector);
};
*/


#endif /* EOS_IG_H_ */
