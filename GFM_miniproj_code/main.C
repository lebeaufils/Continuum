#include "GhostFluid.h"

int main(void){

	//use a switch function here to choose between idealgas or stiffened gas.
	//do the same for choice of solver.
	//switch can be activated through user input or mapping through filename


	gfmTests Tests(1000, 1.0); //(N, L)
	//Tests.testB_Wang();
	//Tests.testB();
	//Tests.testMach10_2();
	//Tests.testMach10_2();
	Tests.switch_test();
	Tests.switch_resolution();

	//Tests.set_EOS(eos1, eos2);
	GhostFluidMethods gfmProblem(0.5, Tests); //See MUSCL.pdf paper forr stability condition suggesting 0.5

	if (Tests.stiffgas == true){
		try{
			gfmProblem.exact_solver_SG(Tests);
		}
		catch (const char* c){
			std::cout << c << std::endl; 
		}
		
		try{
			gfmProblem.solver_RP_SG(Tests);
		}
		catch (const char* c){
			std::cout << c << std::endl; 
		}
	}

	else {
		try{
			gfmProblem.exact_solver(Tests);
		}
		catch (const char* c){
			std::cout << c << std::endl; 
		}
		
		try{
			gfmProblem.solver_RP(Tests);
		}
		catch (const char* c){
			std::cout << c << std::endl; 
		}
	}

}

