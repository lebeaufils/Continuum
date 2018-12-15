/*
 * main.C
 *
 *  Created on: 15 Dec 2018
 *      Author: forte
 */

#include "RPSolvers/HLLC.h"

int main(void){

	//IdealGas EOS;
	JWL EOS;

	eulerTests Tests(100, 1.0);
	Tests.test2();

	HLLC var(0.9, Tests);

	var.initial_conditions(Tests, EOS);
	//var.solver_IG(EOS, Tests);
	var.solver(EOS, Tests);
	var.output(EOS);
}

