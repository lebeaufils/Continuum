#include "headerFiles/Solvers.h"



int main(void){

	eulerTests2D Tests(11);
	Tests.test1();
	//Tests.switch_resolution();
	//Tests.switch_test();

	//EXACT theExact(Tests);
	//theExact.solver(Tests);

	MUSCL::initial_conditions(Tests);
	//MUSCL::solver(Tests.var, 0.9);
	//MUSCL::output(Tests.var);

}

