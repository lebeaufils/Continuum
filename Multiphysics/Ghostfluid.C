#include "GhostFluid.h"

GhostFluidMethods::GhostFluidMethods(double c, gfmTests Test)
	: LevelSetFunction(Test), var1(NULL), var2(NULL), CFL(c), Smax(0), dt(0), count(0){
}

void GhostFluidMethods::initialise_HLLC(gfmTests Test){
	var1 = new HLLC(Test);
	var2 = new HLLC(Test);
}

void GhostFluidMethods::ghost_boundary(RPsolvers* Ureal, EOS* eosReal, RPsolvers* Ughost, EOS* eosGhost, int i){

	//original GFM - constant entropy extrapolation
	//			   continuous Pressure and Velocity

	double dreal = Ureal->U(i, 0); //density
	double velocityreal = Ureal->U(i, 1)/Ureal->U(i, 0);
	double Preal = eosReal->Pressure(Ureal->U, i);

	double Pghost = Preal; 
	double velocityghost = velocityreal;
	double dghost = pow(pow(dreal, eosReal->y)*(Pghost/Preal), 1./eosGhost->y);

	//setting the conservative form of the ghost variables
	Ughost->U(i, 0) = dghost;
	Ughost->U(i, 1) = velocityghost*dghost;
	Ughost->U(i, 2) = Pghost/(eosGhost->y-1) + 0.5*dghost*pow(velocityghost, 2.0);
}

void GhostFluidMethods::ghost_boundary_RP(){
	/*
	double yreal, yghost;
	double Preal, Pghost;
	double dghost;
	double dreal = Ureal(i, 0);
	double velocityreal = Ureal(i, 1)/Ureal(i, 0);
	double velocityghost;

	if (realmaterial == 1){
		yreal = gamma1;
		yghost = gamma2;
	}

	else {
		yreal = gamma2;
		yghost = gamma1;
	}

	Preal = (yreal-1)*(Ureal(i, 2) - 0.5*Ureal(i, 0)*pow((Ureal(i, 1)/Ureal(i, 0)),2.0)); //extrapolated Pressure
	Pghost = Preal;
	velocityghost = velocityreal;
	dghost = pow(pow(dreal, yreal)*(Pghost/Preal), 1.0/yghost);


	//std::cout << Pghost << '\t' << dghost << '\t' << velocityghost << std::endl;

	//setting the conservative form of the ghost variables
	Ughost(i, 0) = dghost;
	Ughost(i, 1) = velocityghost*dghost;
	Ughost(i, 2) = Pghost/(yghost-1) + 0.5*dghost*velocityghost*velocityghost;
	*/
}

//Initial conditions for 2 material system
void GhostFluidMethods::initial_conditions_HLLC(EOS* eos1, EOS* eos2, gfmTests Test){
	initialise_HLLC(Test);

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
			var1->U.row(i+1) = euler1; //setting real values of material 1
		}
		else{
			var2->U.row(i+1) = euler2; //setting real values of material 2
		}
	}
	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=1; i<N+1; i++){
		if (get_sgn(phi(i)) < 0){
			ghost_boundary(var1, eos1, var2, eos2, i); //matrix 2 contains ghost points
		}

		else {
			ghost_boundary(var2, eos2, var1, eos1, i); //matrix 1 contains ghost points
		}
	}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::update_levelset(double dt){
	matrix phi_1; phi_1.resize(N+2, 1); //phi at n+1th timestep
	double velocity;

	for (int i=1; i<N+1; i++){
		if (get_sgn(phi(i)) < 0){ //taking >= 0 as the other material
			velocity = var1->U(i, 1)/var1->U(i, 0); //"advection" velocity for the levelset function, phi
		}
		else {
			velocity = var2->U(i, 1)/var2->U(i, 0);
		}
		phi_1(i) = HJ_FirstOrder(velocity, dt, i);
	}

	phi = phi_1;
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::solver(EOS* eos1, EOS* eos2, gfmTests Test){

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			if (get_sgn(phi(i)) < 0){
				ghost_boundary(var1, eos1, var2, eos2, i); //matrix 2 contains ghost points
			}

			else {
				ghost_boundary(var2, eos2, var1, eos1, i); //matrix 1 contains ghost points
			}
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//compute fluxes at current timestep
		for (int i=0; i<N+1; i++){
			var1->compute_fluxes(eos1, i);
			var2->compute_fluxes(eos2, i);
		}

		//set timestep following CFL conditions with max wavespeed between both materials
		Smax = fmax(var1->Smax, var2->Smax);
		dt = CFL*(dx/Smax); 
		if (t + dt > Test.tstop) dt = Test.tstop - t;

		//evolving the levelset equation
		update_levelset(dt);

		//updating both materials, looping through real and ghost points
		for (int i=1; i<N+1; i++){
			var1->conservative_update_formula(dt, dx, i);
			var2->conservative_update_formula(dt, dx, i);
		}

		var1->boundary_conditions();
		var2->boundary_conditions();
		
		t += dt;
		count += 1;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
}

void GhostFluidMethods::output(EOS* eos1, EOS* eos2){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	double d, u, P, e;
	for (int i=1; i<N+1; i++){
		if (get_sgn(phi(i)) < 0){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos1->Pressure(var1->U, i);
			e = eos1->internalE(var1->U, i);
		}

		else {
			d = var2->U(i, 0);
			u = var2->U(i, 1)/var2->U(i, 0);
			P = eos2->Pressure(var2->U, i);
			e = eos2->internalE(var2->U, i);
		}
		outfile << X(i) << '\t' << d << '\t' << u
		<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done: GFM-HLLC" << std::endl;
}





