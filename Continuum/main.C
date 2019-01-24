/*
 * main.C
 *
 *  Created on: 15 Dec 2018
 *      Author: forte
 */

#include "RPSolvers/Solvers.h"

int main(void){

	/*
	IdealGas EOS;
	//JWL EOS;

	eulerTests Tests(100, 1.0); //(N, L)
	Tests.test5();

	FORCE var(0.9, Tests);

	var.initial_conditions(EOS, Tests);

	var.solver(EOS, Tests);
	var.output(EOS);*/

	IdealGas EOS;
	eulerTests Tests(100, 1.0); //(N, L)
	Tests.test1();

	RPsolvers *var = new MUSCL(0.9, Tests);

	//HLLCeuler.boundary_conditions();
	var->initial_conditions(EOS, Tests);
	var->solver(EOS, Tests);
	var->output(EOS);
}

