#include "FVschemes.H"

#include <fstream>
#include <iomanip>

//constructor
FiniteVolume::FiniteVolume(int n, double l, double a, double c)
	: N(n), L(l), a(a), c(c), X(N+1, 0), dx(l/n), dt(c*l/(n*a)){
}
scalar::scalar(int n, double l, double a, double c)
	: FiniteVolume(n, l, a, c), u(N+2, 0), u_1(N+2, 0), f(N+1,0){
}

scalarFORCE::scalarFORCE(int n, double l, double a, double c)
	: FiniteVolume(n, l, a, c), scalar(n, l, a, c){
}

linearadvection::linearadvection(int n, double l, double a, double c)
	: FiniteVolume(n, l, a, c), scalar(n, l, a, c), scalarFORCE(n, l, a, c){
}

euler::euler(int n, double l, double a, double c)
	: FiniteVolume(n, l, a, c), U(N+2, std::vector<double>(3)), U_1(N+2, std::vector<double>(3)),
	  F(N+1, std::vector<double>(3)){
}

eulerFORCE::eulerFORCE(int n, double l, double a, double c)
	: FiniteVolume(n, l, a, c), euler(n, l, a, c){
}

//scalar
void scalar::initial_conditions_step()
{
	for (int i=0; i<(N+1); i++){ //fixed boundary. check finite D
		X[i] = i*dx;

		if (X[i] < (L/2)){u[i+1] = 1.0;} //u starts from 1 to N+1 (index N)
		else {u[i] = 0.0;}
	}
	boundary_conditions();
}

void scalar::initial_conditions_square()
{
	for (int i=0; i<(N+1); i++){ //fixed boundary. check finite D
		X[i] = i*dx;

		if (X[i] < 0.3){u[i+1] = 0.0;} //u starts from 1 to N+1 (index N)
		else if (X[i] > 0.7){u[i+1] = 0.0;}
		else {u[i+1] = 1.0;}
	}
	boundary_conditions();
}


void scalar::boundary_conditions()
{
	u[0] = u[1];
	u[N+1] = u[N+2];
}

void scalar::output(std::string filename)
{
	std::ofstream outfile;
	outfile.open(filename);

	for (int i=1; i<(N+2); i++){
		outfile << X[i] << '\t' << u[i] << std::endl;
	}
}

void scalar::plotname(double t)
{
	std::stringstream filename;
	filename << "Plot_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}

//First OrdeR CEntered scheme (scalar)
	// combines 2 schemes
void scalarFORCE::solver()
{
	for (int i=1; i<(N+2); i++){
		u_1[i] = u[i] + (dt/dx)*(f[i-1]-f[i]);
	}
	for (int i=1; i<N+2; i++){
		u[i] = u_1[i];
	}
}

double linearadvection::flux(double q)
{
	return a*q;
}

void linearadvection::evaluate_flux()
{
//LxF
	double LxF;
	double R; //Richtmyer

	for (int i=0; i<(N+1); i++)
	{
		LxF = ((1+c)/(2*c))*flux(u[i]) + ((c-1)/(2*c))*flux(u[i+1]);
		R = flux((0.5*(u[i] + u[i+1]) + 0.5*(dt/dx)*(flux(u[i]) - flux(u[i+1])))); //linadv
		f[i] = 0.5*(LxF + R); //f[i-0.5]
	}
}

void linearadvection::plotname(double t)
{
	std::stringstream filename;
	filename << "linadvFORCE_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}

//Euler
void euler::initial_conditions_step()
{
	double initialL[] = {2,0,10}; //density, speed, pressure
	double initialR[] = {1,0,1};

	for (int i=0; i<(N+1); i++){ //fixed boundary. check finite D
		X[i] = i*dx;

		if (X[i] < (L/2)){
			U[i+1].assign(initialL, &initialL[3]);
		}
		 //u starts from 1 to N+1 (index N)
		else {
			U[i+1].assign(initialR, &initialR[3]);
		}
	}
	boundary_conditions();
}


void euler::boundary_conditions()
{
	U[0] = U[1];
	U[N+1] = U[N+2];
}

void euler::output(std::string filename)
{
	std::ofstream outfile;
	outfile.open(filename);

	for (int i=1; i<(N+2); i++){
		outfile << X[i] << '\t' << U[i][0] << '\t' << U[i][1] << '\t' << U[i][2] << std::endl;
	}
}

void euler::plotname(double t)
{
	std::stringstream filename;
	filename << "Euler_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}


