#include "GhostFluid.h"

int main(void){

	/*//use a switch function here to choose between idealgas or stiffened gas.
	//do the same for choice of solver.
	//switch can be activated through user input or mapping through filename

 	EOS* eos1 = new IdealGas();
	EOS* eos2 = new IdealGas();
	EOS* eos3 = new IdealGas();

	gfmTests Tests(100, 1.0); //(N, L)

	Tests.testB();
	//Tests.test_example_1();

	GhostFluidMethods gfmProblem(0.5, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5
	//gfmProblem.initial_conditions_HLLC(eos1, eos2, Tests);

	try{
		//gfmProblem.initial_conditions_RP(eos1, eos2, eos3, Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	//gfmProblem.solver_RP(Tests);
	try{
		gfmProblem.exact_solver(Tests, eos1, eos2, eos3);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	//gfmProblem.output(eos1, eos2, eos3);

	delete eos1; delete eos2; delete eos3;

*/
	//Note no exact solver exists for such a configuration

/*	gfmTests Tests(400, 1.0); //(N, L)

	//Tests.testB_Wang();
	Tests.testSG();

	GhostFluidMethods gfmProblem(0.5, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5

	try{
		gfmProblem.solver_RP(Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	try{
		gfmProblem.exact_solver(Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
*/	


	gfmTests Tests(100, 1.0); //(N, L)
	//Tests.testB_Wang();
	//Tests.testB();
	Tests.testMach10();

	//Tests.set_EOS(eos1, eos2);

	GhostFluidMethods gfmProblem(0.5, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5

	try{
		gfmProblem.exact_solver(Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}
	
	try{
		gfmProblem.solver(Tests);
	}
	catch (const char* c){
		std::cout << c << std::endl; 
	}

}

