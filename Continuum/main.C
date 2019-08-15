#include "headerFiles/Solvers.h"
#include "headerFiles/LevelSet.h"
#include "headerFiles/RigidBodies.h"
#include <typeinfo>


int main(void){


	demTests Tests(401);
	Tests.test1();

	RigidBodies::rigid_body_solver(Tests, 0.5);
}

