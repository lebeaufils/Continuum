#include "headerFiles/Solvers.h"



int main(void){

	eulerTests2D Tests(101);
	Tests.test3();


	//eulerTests Test1d(101);
	//Test1d.test1_stationary();
	//Tests.switch_resolution();
	//Tests.switch_test();

	//EXACT theExact(Test1d);
	//theExact.solver(Test1d);

	MUSCL::initial_conditions(Tests);
	MUSCL::solver(Tests.var, 0.2);
	MUSCL::output(Tests.var);

}

