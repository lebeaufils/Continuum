/*
 * main.C
 *
 *  Created on: 15 Dec 2018
 *      Author: forte
 */

#include "RPSolvers/MUSCL.h"

int main(void){

	IdealGas EOS;
	//JWL EOS;

	eulerTests Tests(100, 1.0);
	Tests.test1();

	MUSCL var(0.9, Tests);

	var.initial_conditions(EOS, Tests);
	//var.solver_IG(EOS, Tests);
	var.solver(EOS, Tests);
	var.output(EOS);
}

