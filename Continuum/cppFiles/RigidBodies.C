#include "../headerFiles/RigidBodies.h"

void RigidBodies::rigid_boundary(Euler1D real, Euler1D &ghost, int i){
//Reflective boundary for a rigid body
	//cell i is the cell just before the rigid body in i+1
	double density = real.U(i, 0); //density and pressure remain constant across the rigid body interface
	double pressure = real.state_function->Pressure(real.U, i);
	double velocity = -real.U(i, 1)/real.U(i, 0); //reverse direction in normal velocity
	
	//Primitive form of the ghost variables
	vector primitive(density, velocity, pressure);

	//Redefining the real fluids at the boundary (isobaric fix, changing only density)
	for (int j=0; j<3; j++){
		//Defining the states at ghost fluid cells [i+1, i+2, i+3],
		//checking that the ghost cells remain within computational domain
		if(i+j+1 < real.N+2) ghost.U.row(i+j+1) = real.state_function->conservedVar(primitive); 
	}
}

void RigidBodies::initial_conditions(LevelSetFunction1D levelset, rigidTests &Test){

}

void RigidBodies::solver(){

}

void RigidBodies::output(){

}

void RigidBodies::rigid_body_solver();