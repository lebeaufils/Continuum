/*----------------------------------------------------------------------------------
 Constant Volume explosion.
 One cell that is filled with reactive gas.
 "Heating" this cell gives the thermal runaway in 0th dimension
 This is a rapid increase in temperature and drop in species concentration
 ----------------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>

//Initial conditions
const double T0 = 1.0;
const double c0 = 1.0;

const double tstop = 2.0;

//constants
const double Q = 1.24;
const double Ea = 14.99;
const double epsilon = 1/Ea;
const double y = 1.4;

void semi_analytic_solver(double &c, double &T, double dt){

	double c_1, T_1;
	//double epsilon = 1/Ea;
	double beta = Q*y*c0/T0;

	c_1 = c*exp(-dt*(c*epsilon*T0/beta)*exp((1/(epsilon*T0))*(1-1/T)));
	T_1 = 1 + beta*(1-c_1);

	//update c and T
	c = c_1;
	T = T_1;
}

double fT(double c, double T){
	return c*epsilon*T0*exp(1./(epsilon*T0)*(1. - 1./T));
}

double fc(double c, double T){
	double beta = Q*y*c0/T0;
	return -c*epsilon*T0/beta*exp(1./(epsilon*T0)*(1. - 1./T));
}

void RK4(double &c, double &T, double dt){

	double Tk1 = dt*fT(c, T);
	double ck1 = dt*fc(c, T);
	double Tk2 = dt*fT(c + ck1/2, T + Tk1/2);
	double ck2 = dt*fc(c + ck1/2, T + Tk1/2);
	double Tk3 = dt*fT(c + ck2/2, T + Tk2/2);
	double ck3 = dt*fc(c + ck2/2, T + Tk2/2);
	double Tk4 = dt*fT(c + ck3, T + Tk3);
	double ck4 = dt*fc(c + ck3, T + Tk3);

	double c_1 = c + 1./6*(ck1 + 2*ck2 + 2*ck3 + ck4);
	double T_1 = T + 1./6*(Tk1 + 2*Tk2 + 2*Tk3 + Tk4);

	c = c_1;
	T = T_1;
}

int main(void){

	/*
	std::ofstream outfile;
	outfile.open("ODEsolver1.txt");

	double c = c0;
	double T = T0;
	double t = 0;
	double dt = 0.001;

	outfile << t << '\t' << c << '\t' << T << std::endl;
	std::cout << t << '\t' << c << '\t' << T << std::endl;

	int count = 0;
	do{
		if (t + dt > tstop) dt = tstop - t;
		t += dt;

		semi_analytic_solver(c, T, dt);
		outfile << t << '\t' << c << '\t' << T << std::endl;
		std::cout << t << '\t' << c << '\t' << T << std::endl;

		count += 1;
	}while (t < tstop);
	std::cout << count << std::endl;
	outfile.close();
	*/

	std::ofstream outfile;
	outfile.open("RK4solver01.txt");

	double c = c0;
	double T = T0;
	double t = 0;
	double dt = 0.01;

	outfile << t << '\t' << c << '\t' << T << std::endl;
	std::cout << t << '\t' << c << '\t' << T << std::endl;

	int count = 0;
	do{
		if (t + dt > tstop) dt = tstop - t;
		t += dt;

		RK4(c, T, dt);
		outfile << t << '\t' << c << '\t' << T << std::endl;
		std::cout << t << '\t' << c << '\t' << T << std::endl;

		count += 1;
	}while (t < tstop);
	std::cout << count << std::endl;
	outfile.close();
}
