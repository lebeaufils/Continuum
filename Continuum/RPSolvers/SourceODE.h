/*
 * SourceODE.h
 *
 *  Created on: 22 Dec 2018
 *      Author: forte
 */

#ifndef SOURCEODE_H_
#define SOURCEODE_H_

typedef double (*myFunc)(double x, double y);

class ODEsolvers
{
	double c, T; //f(c, T)

public:
	ODEsolvers(double, double); //initial values x0 and y0

	void RK4(double, myFunc, myFunc);
	void RK4(double, double&, double&, myFunc);
	//void semi_analytic(double);
	//void output();
};




#endif /* SOURCEODE_H_ */
