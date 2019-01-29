#include "GhostFluid.h"

GhostFluidMethods::GhostFluidMethods(gfmTests Test)
	: LevelSetFunction(Test), Smax(0), dt(0){
	//Using "T" solver
	var1 = NULL;
	var2 = NULL;
}


void GhostFluidMethods::ghost_boundary(gfmTests Test, RPsolvers* Ureal, EOS* eosReal, RPsolvers* Ughost, EOS* eosGhost, int i){
/*
	original GFM - constant entropy extrapolation
				   continuous Pressure and Velocity
*/
	double dreal = Ureal->U(i, 0); //density
	double velocityreal = Ureal->U(i, 1)/Ureal->U(i, 0);
	double Preal = eosReal.Pressure(U, i);

	double Pghost = Preal; 
	double velocityghost = velocityreal;
	double dghost = pow(pow(dreal, eosReal->y)*(Pghost/Preal), 1./eosGhost->y);

	//setting the conservative form of the ghost variables
	Ughost->U(i, 0) = dghost;
	Ughost->U(i, 1) = velocityghost*dghost;
	Ughost->U(i, 2) = Pghost/(eosGhost.y-1) + 0.5*dghost*pow(velocityghost, 2.0);
}

void GhostFluidMethods::ghost_boundary_RP(){

}

//Initial conditions for 2 material system
void GhostFluidMethods::initial_conditions<HLLC>(HLLC varL, HLLC varR, EOS* eos1, EOS* eos2, gfmTests Test){
	var1 = &varL;
	var2 = &varR;

	vector euler1;
	vector euler2;

	//setting material EOS parameter
	eos1->y = Test.y1;
	eos2->y = Test.y2;

	//initial values for conserved variables
	euler1 = eos1->conservedVar(Test.initialL);
	euler2 = eos2->conservedVar(Test.initialR);

	for (int i=0; i<N; i++){
		X(i+1) = i*dx;
		signed_distance_function_1D(i);
		if (X(i+1)  < x0){
			var1->row(i+1) = euler1; //setting real values of material 1
		}
		else{
			var2->row(i+1) = euler2; //setting real values of material 2
		}
	}
	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=0; i<N+1; i++){
		if (X(i+1) < x0){
			ghost_boundary(Test, var1, eos1, var2, eos2, i) //matrix 2 contains ghost points
		}

		else {
			ghost_boundary(Test, var2, eos2, var1, eos1, i) //matrix 1 contains ghost points
		}
	}

	var1->boundary_conditiona();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::update_level_set(double dt){
	matrix phi_1; phi_1.resize(N+2, 1); //phi at n+1th timestep
	double velocity;

	for (int i=1; i<N+1; i++){
		if (get_sgn(phi(i)) < 0){ //taking >= 0 as the other material
			velocity = var1->U(i, 1)/var1->U(i, 0); //"advection" velocity for the levelset function, phi
		}
		else {
			velocity = var2->U(i, 1)/var2->U(i, 0);
		}
		phi_i(i) = HJ_FirstOrder(velocity, dt, i);
	}

	phi = phi_1;
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::solver(gfmTests test){

	
	double t = 0.0;
	do{

	//compute fluxes at current timestep
	for (int i=0; i<N+1; i++){
		var1->compute_fluxes_HLLC(eos1, i);
		var2->compute_fluxes_HLLC(eos2, i);
	}

	update_level_set(dt);

	//set timestep following CFL conditions with max wavespeed between both materials
	Smax = max(var1->Smax, var2->Smax);
	dt = CFL*(dx/Smax); 
	if (t + dt > Test.tstop) dt = Test.tstop - t;

	//updating var1 and var2
	for (int i=1; i<N+1; i++){
		var1->conservative_update_formula();
	}
	

	}while(t<test.tstop)






}




