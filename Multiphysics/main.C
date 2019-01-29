/*
 * main.C
 *
 *  Created on: 15 Dec 2018
 *      Author: forte
 */

#include "RPsolvers/Solvers.h"

int main(void){

	//use a switch function here to choose between idealgas or stiffened gas.
	//do the same for choice of solver.
	//switch can be activated through user input or mapping through filename
	EOS *SG = new IdealGas();
	//SG->GetGamma();

	eulerTests Tests(100, 1.0); //(N, L)
	Tests.test1();

	RPsolvers *var = new HLLC(0.9, Tests);

	//HLLCeuler.boundary_conditions();
	var->initial_conditions(SG, Tests);
	var->solver(SG, Tests);
	var->output(SG);

	delete SG;
	delete var;
}

