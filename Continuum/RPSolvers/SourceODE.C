/*
 * SourceODE.C
 *
 *  Created on: 22 Dec 2018
 *      Author: forte
 */
#include "SourceODE.h"

ODEsolvers::ODEsolvers(double c0, double T0)
	: c(c0), T(T0) {}

void ODEsolvers::RK4(double h, double &x, double &y, double (*f)(double, double)){ //f = f(x, y)
	double k1 = h*f(x, y);
	double k2 = h*f(x + h/2, y + k1/2);
	double k3 = h*f(x + h/2, y + k2/2);
	double k4 = h*f(x + h, y + k3);

	double y_1 = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4); // + O(h^5)
}

void ODEsolvers::RK4(double dt, double (*fc)(double, double), double (*fT)(double, double)){
	//double (*fc)(double c, double T);
	//fc = &my_function;
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


