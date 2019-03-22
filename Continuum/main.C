#include "RPsolvers/Solvers.h"



int main(void){

	eulerTests Tests(100);
	Tests.test1();
	Tests.switch_resolution();
	Tests.switch_test();

	EXACT theExact(Tests);
	theExact.solver(Tests);

	MUSCL::initial_conditions_1D(Tests);
	MUSCL::solver(Tests.var, 0.9);
	MUSCL::output(Tests.var);

}

