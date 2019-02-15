#include "GhostFluid.h"

int main(void){

	//use a switch function here to choose between idealgas or stiffened gas.
	//do the same for choice of solver.
	//switch can be activated through user input or mapping through filename

 /*	EOS* eos1 = new IdealGas();
	EOS* eos2 = new IdealGas();
	EOS* eos3 = new IdealGas();

	gfmTests Tests(400, 1.0); //(N, L)

	Tests.testB();
	//Tests.test_example_1();

	GhostFluidMethods gfmProblem(0.5, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5
	//gfmProblem.initial_conditions_HLLC(eos1, eos2, Tests);

	try{
		gfmProblem.initial_conditions(eos1, eos2, eos3, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	gfmProblem.solver(eos1, eos2, eos3, Tests);
	try{
		gfmProblem.exact_solver(Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	gfmProblem.output(eos1, eos2, eos3);

	delete eos1; delete eos2; delete eos3;
*/

/*	//Note no exact solver exists for such a configuration
	EOS* eos1 = new IdealGas();
	EOS* eos2 = new IdealGas();
	EOS* eos3 = new IdealGas();
	EOS* eos4 = new IdealGas();

	gfmTests Tests(400, 1.0); //(N, L)

	Tests.testB_Wang();

	GhostFluidMethods gfmProblem(0.6, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5
	//gfmProblem.initial_conditions_HLLC(eos1, eos2, Tests);
	try{
		gfmProblem.initial_conditions(eos1, eos2, eos3, eos4, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	gfmProblem.solver(eos1, eos2, eos3, eos4, Tests);
	gfmProblem.output(eos1, eos2, eos3, eos4);

	delete eos1; delete eos2; delete eos3; delete eos4;
*/
	//EOS* eos1 = new IdealGas();
	//EOS* eos2 = new IdealGas();

	EOS* eos1 = new StiffenedGas();
	EOS* eos2 = new StiffenedGas();

	EOS* eos3 = new IdealGas();
	EOS* eos4 = new IdealGas();

	gfmTests Tests(400, 1.0); //(N, L)
	Tests.testF();
	//Tests.test_example_1();

	Tests.set_EOS(eos1, eos2);
	Tests.set_EOS(eos3, eos4);

	vector WL = Tests.initialL;	
	vector WR = Tests.initialR;


	eos1->y_constants(WL);
	eos2->y_constants(WR);

	eos3->y_constants(WL);
	eos4->y_constants(WR);


	GhostFluidMethods gfmProblem(0.5, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5
	
	std::cout << gfmProblem.compute_star_pressure(eos3, eos4) << std::endl;
	std::cout << gfmProblem.compute_star_pressure_SG(eos1, eos2) << std::endl;

	try{
		gfmProblem.exact_solver_SG(Tests, eos1, eos2);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	
	/*try{
		//gfmProblem.initial_conditions_RP(eos1, eos2, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	try{
		//gfmProblem.solver_RP(eos1, eos2, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	try{
		//gfmProblem.initial_conditions(eos1, eos2, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}

	try{
		//gfmProblem.solver(eos1, eos2, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	gfmProblem.output(eos1, eos2);
	try{
		gfmProblem.exact_solver(Tests, eos1, eos2);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	*/
	delete eos1; delete eos2; delete eos3; delete eos4;
/*
	//note exact solver for stiffened gas not written yet
	EOS* IG = new IdealGas();
	eulerTests Test(400, 1.0);
	//Test.test1_stationary();
	Test.test6();

	EXACT var(Test, IG);
	var.initial_conditions();

	var.sampling(Test.tstop);
	var.output();

	delete IG;

    return 0;
 */
}

