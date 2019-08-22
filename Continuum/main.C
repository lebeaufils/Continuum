#include "headerFiles/Solvers.h"
#include "headerFiles/LevelSet.h"
#include "headerFiles/RigidBodies.h"
#include <typeinfo>


int main(void){

	//Generate the domain with the demTests constructor.
	demTests Tests(801,401);
	Tests.test6(); //initialise the test case (1-6)
	//Call the solver using parameters (Test, CFL)
	RigidBodies::rigid_body_solver(Tests, 0.5);
}

