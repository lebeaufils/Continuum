/*
 * main.C
 *
 *  Created on: 15 Dec 2018
 *      Author: forte
 */

#include "GhostFluid.h"

int main(void){

	//use a switch function here to choose between idealgas or stiffened gas.
	//do the same for choice of solver.
	//switch can be activated through user input or mapping through filename
	/*
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
	*/

/*
	EOS* eos1 = new IdealGas();
	EOS* eos2 = new IdealGas();

	gfmTests Tests(400, 1.0); //(N, L)

	Tests.testA();
	//Tests.test_example_1();

	GhostFluidMethods gfmProblem(0.9, Tests);
	gfmProblem.initial_conditions_HLLC(eos1, eos2, Tests);
	gfmProblem.solver(eos1, eos2, Tests);
	gfmProblem.output(eos1, eos2);

	delete eos1; delete eos2;
*/
	EOS* IG = new IdealGas();
	eulerTests Test(100, 1.0);
	Test.test4();

	EXACT var(Test);
	var.initial_conditions(Test);
	double pstar = var.compute_star_pressure(Test.initialL, Test.initialR, IG);
	std::cout << pstar << std::endl;

	delete IG;
}

