#include "headerFiles/Solvers.h"



int main(void){

	//eulerTests2D Tests(101);
	//Tests.test3();
	//MUSCL::muscl_solver(Tests, 0.9);

	//eulerTests Test1d(101);
	//Test1d.test1_stationary();
	//Tests.switch_resolution();
	//Tests.switch_test();

	//EXACT theExact(Test1d);
	//theExact.solver(Test1d);

	rigidTests Tests(11);
	Tests.test1();
}

