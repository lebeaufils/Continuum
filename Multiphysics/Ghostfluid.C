#include "GhostFluid.h"

GhostFluidMethods::GhostFluidMethods(double c, gfmTests Test)
	: LevelSetFunction(Test), var1(NULL), var2(NULL), CFL(c), Smax(0), dt(0), count(0){
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

	//std::cout << dreal << '\t' << velocityreal << '\t' << Preal << std::endl;
}

void GhostFluidMethods::ghost_boundary_RP(){
//Exact riemannsolver for idealgas and stiffened gas. Ghost point values are the RP star states.
}

//--------------------------------------------------------
//	HLLC
//--------------------------------------------------------

//Initial conditions for 2 material system
void GhostFluidMethods::initial_conditions_HLLC(EOS* eos1, EOS* eos2, gfmTests Test){
	var1 = new HLLC(Test);
	var2 = new HLLC(Test);

	vector euler1;
	vector euler2;

	//setting material EOS parameter
	eos1->y = Test.y1;
	eos2->y = Test.y2;

	//initial values for conserved variables
	euler1 = eos1->conservedVar(Test.initialL);
	euler2 = eos2->conservedVar(Test.initialR);


	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
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

//Initial conditions for 3 material system
void GhostFluidMethods::initial_conditions_HLLC(EOS* eos1, EOS* eos2, EOS* eos3, gfmTests Test){
	var1 = new HLLC(Test);
	var2 = new HLLC(Test);

	//setting material EOS parameter
	eos1->y = Test.y1;
	eos2->y = Test.y2;
	if (Test.y3 != 0){
		eos3->y = Test.y3;
	}
	else {
		throw "Number of materials do not match the chosen solver.";
	}


	//initial values for conserved variables
	vector euler1 = eos1->conservedVar(Test.initialL);
	vector euler2 = eos2->conservedVar(Test.initialM1);
	vector euler3 = eos3->conservedVar(Test.initialR);

	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		signed_distance_function_1D_2(i);
		if (X(i+1)  <= x0){
			var1->U.row(i+1) = euler1; //setting real values of Left initial condition
		}
		else if (X(i+1) > x1){
			var1->U.row(i+1) = euler3; //setting real values of Right initial condition
		}
		else{
			var2->U.row(i+1) = euler2; //Middle initial condition
		}
	}

	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=1; i<N+1; i++){
		int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
		//std::cout << phi(i) << '\t' << testsgn << std::endl;
		if (testsgn > -2 && testsgn < 1 && i!=N){ //-1 or 0, the boundary is crossed
			if (get_sgn(phi(i)) < 0){//the interface is to the right of cell i
				for (int j=0; j<3; j++){
					//to the left of interface, the ghostfluid is in var2
					ghost_boundary(var1, eos1, var2, eos2, i-j); 
					//to the right of the interface, the ghostfluid is in var1
					ghost_boundary(var2, eos2, var1, eos1, i+j+1);
					//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
					//i, i-1 and i-2, while thst of var1 is i+1, i+2 and i+3.
				}
			}
			else {//the interface is to the left of cell i+1 (i+1 is negative while i is positive)
				for(int j=0; j<3; j++){
					//to the left of interface, the ghostfluid is in var1
						//note that we are now in the right material which has EOS: eos3
					ghost_boundary(var2, eos2, var1, eos3, i-j); 
					//to the right of the interface, the ghostfluid is in var2
					ghost_boundary(var1, eos3, var2, eos2, i+1+j);
					//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
					//i, i+1 and i+2, while thst of var1 is i-1, i-2 and i-3.
				}
			}
		} 
	}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::update_levelset_HLLC(double dt){
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

void GhostFluidMethods::solver_HLLC(EOS* eos1, EOS* eos2, gfmTests Test){

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			if (get_sgn(phi(i-1)) < 0){
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
		update_levelset_HLLC(dt);

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

void GhostFluidMethods::solver_HLLC(EOS* eos1, EOS* eos2, EOS* eos3, gfmTests Test){

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i<N-2){ //-1 or 0, the boundary is crossed
				if (get_sgn(phi(i)) < 0){//the interface is to the right of cell i
					for (int j=0; j<3; j++){
						//to the left of interface, the ghostfluid is in var2
						ghost_boundary(var1, eos1, var2, eos2, i-j); 
						//to the right of the interface, the ghostfluid is in var1
						ghost_boundary(var2, eos2, var1, eos1, i+j+1);
						//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
						//i, i-1 and i-2, while thst of var1 is i+1, i+2 and i+3.
					}
				}
				else {//the interface is to the left of cell i+1 (i+1 is negative while i is positive)
					for(int j=0; j<3; j++){
						//to the left of interface, the ghostfluid is in var1
							//note that we are now in the right material which has EOS: eos3
						ghost_boundary(var2, eos2, var1, eos3, i-j); 
						//to the right of the interface, the ghostfluid is in var2
						ghost_boundary(var1, eos3, var2, eos2, i+1+j);
						//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
						//i, i+1 and i+2, while thst of var1 is i-1, i-2 and i-3.
					}
				}
			} 

			if (testsgn > -2 && testsgn < 1 && i>N-2){
				std::cout << "Warning, interface is going out of bounds" << std::endl;
			}
			//std::cout << var1->U.row(i) << '\t' << var2->U.row(i) << std::endl;
		}

		var1->boundary_conditions();
		var2->boundary_conditions();
		/*for (int i=0; i<N+2; i++){
			if (count ==0 ) std::cout << var1->U.row(i) << std::endl;
		}*/

		//compute fluxes at current timestep
		
		/*for (int i=0; i<N+1; i++){
			var1->compute_fluxes(eos1, i);
			var2->compute_fluxes(eos2, i);
		}*/

		//compute fluxes at current timestep
		
		bool flag = true;
		//std::cout << "test1" << std::endl;
		for (int i=0; i<N+1; i++){
			if (phi(i) < 0 && flag==true){
				//std::cout <<"a " << i << std::endl;
				var1->compute_fluxes(eos1, i);
				//if (get_sgn(phi(i+1)) > 0) var1->compute_fluxes(eos1, i+1);
			}
			if (phi(i) >= 0){
				//std::cout <<"b " << i << std::endl;
				if (get_sgn(phi(i-1)) < 0) var2->compute_fluxes(eos1, i-1);
				var2->compute_fluxes(eos2, i);
				//if (get_sgn(phi(i+1)) < 0) var1->compute_fluxes(eos1, i+1);
				flag = false;
			}
			if (phi(i) < 0 && flag == false){
				//std::cout <<"c " << i << std::endl;
				if (get_sgn(phi(i-1)) > 0) var1->compute_fluxes(eos3, i-1);
				var1->compute_fluxes(eos3, i);
			}
			//if (count == 0 ) std::cout << var1->U.row(i) << std::endl; //'\t' << var2->U.row(i) << std::endl;
			//if (count ==1 )std::cout << i << '\t' << var1->F.row(i) << std::endl;//'\t' << var2->F.row(i) << std::endl;
		//note error of computing fluxes when the boundary lies exacdtly at cell i!!
		}
		
		//set timestep following CFL conditions with max wavespeed between both materials
		Smax = fmax(var1->Smax, var2->Smax);
		dt = CFL*(dx/Smax); 
		if (t + dt > Test.tstop) dt = Test.tstop - t;



		//updating both materials, looping through real and ghost points
		for (int i=1; i<N+1; i++){
			if (phi(i) < 0){
				//if (get_sgn(phi(i+1)) > 0) var1->conservative_update_formula(dt, dx, i+1);
				//if (get_sgn(phi(i-1)) > 0) var1->conservative_update_formula(dt, dx, i-1);
				//if (count ==0 )std::cout << i << '\t' << var1->F.row(i-1) << std::endl;
				var1->conservative_update_formula(dt, dx, i);
			}
			else {
				//if (get_sgn(phi(i-1)) < 0) var2->conservative_update_formula(dt, dx, i-1);
				var2->conservative_update_formula(dt, dx, i);
				//if (get_sgn(phi(i+1)) < 0) var2->conservative_update_formula(dt, dx, i+1);
			}
			//if (phi(i)<0) var1->conservative_update_formula(dt, dx, i);
			//else var2->conservative_update_formula(dt, dx, i);
		}

		//evolving the levelset equation
		update_levelset_HLLC(dt);

		var1->boundary_conditions();
		var2->boundary_conditions();
		
		t += dt;
		count += 1;
		//if (count%50 ==0) std::cout << count << std::endl;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
}

void GhostFluidMethods::output_HLLC(EOS* eos1, EOS* eos2){
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

void GhostFluidMethods::output_HLLC(EOS* eos1, EOS* eos2, EOS* eos3){
	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	double d, u, P, e;
	bool flag = true;
	for (int i=1; i<N+1; i++){
		if (phi(i) < 0 && flag == true){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos1->Pressure(var1->U, i);
			e = eos1->internalE(var1->U, i);
		}

		if (phi(i) >= 0) {
			d = var2->U(i, 0);
			u = var2->U(i, 1)/var2->U(i, 0);
			P = eos2->Pressure(var2->U, i);
			e = eos2->internalE(var2->U, i);
			flag = false;
		}
		if (phi(i) < 0 && flag == false){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos3->Pressure(var1->U, i);
			e = eos3->internalE(var1->U, i);
		}

		outfile << X(i) << '\t' << d << '\t' << u
		<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done: GFM-HLLC" << std::endl;
}

//--------------------------------------------------------
//	MUSCL
//--------------------------------------------------------

void GhostFluidMethods::initial_conditions_MUSCL(EOS* eos1, EOS* eos2, gfmTests Test){
	var1 = new MUSCL(Test);
	var2 = new MUSCL(Test);

	vector euler1;
	vector euler2;

	//setting material EOS parameter
	eos1->y = Test.y1;
	eos2->y = Test.y2;

	//initial values for conserved variables
	euler1 = eos1->conservedVar(Test.initialL);
	euler2 = eos2->conservedVar(Test.initialR);


	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		signed_distance_function_1D(i); //refers to x(i+1)
		if (X(i+1)  <= x0){
			var1->U.row(i+2) = euler1; //setting real values of material 1
		}
		else{
			var2->U.row(i+2) = euler2; //setting real values of material 2
		}
	}
	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=0; i<N; i++){
		if (get_sgn(phi(i+1)) < 0){
			ghost_boundary(var1, eos1, var2, eos2, i+2); //matrix 2 contains ghost points
		}

		else {
			ghost_boundary(var2, eos2, var1, eos1, i+2); //matrix 1 contains ghost points
		}
	}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction

}

void GhostFluidMethods::initial_conditions_MUSCL(EOS* eos1, EOS* eos2, EOS* eos3, gfmTests Test){
	var1 = new MUSCL(Test);
	var2 = new MUSCL(Test);

		//setting material EOS parameter
	eos1->y = Test.y1;
	eos2->y = Test.y2;
	if (Test.y3 != 0){
		eos3->y = Test.y3;
	}
	else {
		throw "Number of materials do not match the chosen solver.";
	}

	//initial values for conserved variables
	vector euler1 = eos1->conservedVar(Test.initialL);
	vector euler2 = eos2->conservedVar(Test.initialM1);
	vector euler3 = eos3->conservedVar(Test.initialR);

	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		signed_distance_function_1D_2(i);
		if (X(i+1)  <= x0){
			var1->U.row(i+2) = euler1; //setting real values of Left initial condition
		}
		else if (X(i+1) > x1){
			var1->U.row(i+2) = euler3; //setting real values of Right initial condition
		}
		else{
			var2->U.row(i+2) = euler2; //Middle initial condition
		}

		//std::cout << var1->U.row(i) << '\t' << var2->U.row(i) << std::endl;
	}

	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=1; i<N+1; i++){
		int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
		//std::cout << phi(i) << '\t' << testsgn << std::endl;
		if (testsgn > -2 && testsgn < 1 && i!=N){ //-1 or 0, the boundary is crossed
			if (get_sgn(phi(i)) < 0){//the interface is to the right of cell i
				for (int j=0; j<4; j++){
					//to the left of interface, the ghostfluid is in var2
					ghost_boundary(var1, eos1, var2, eos2, (i+1)-j); 
					//to the right of the interface, the ghostfluid is in var1
					ghost_boundary(var2, eos2, var1, eos1, (i+1)+j+1);
					//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
					//i, i-1 and i-2, while thst of var1 is i+1, i+2 and i+3.
				}
			}
			else {//the interface is to the left of cell i+1 (i+1 is negative while i is positive)
				for(int j=0; j<4; j++){
					//std::cout << var2->U.row(i+1-j) << std::endl;
					//to the left of interface, the ghostfluid is in var1
						//note that we are now in the right material which has EOS: eos3
					ghost_boundary(var2, eos2, var1, eos3, (i+1)-j);
					//to the right of the interface, the ghostfluid is in var2
					ghost_boundary(var1, eos3, var2, eos2, (i+1)+1+j);
					//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
					//i, i+1 and i+2, while thst of var1 is i-1, i-2 and i-3.
				}
			}
		} 
	}

	/*`for (int i=1; i<N+1; i++){
		std::cout << i << '\t' << get_sgn(phi(i)) << '\t' <<var1->U.row(i+1) << '\t' << var2->U.row(i+1) << std::endl;
	}*/

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction

}

void GhostFluidMethods::update_levelset_MUSCL(double dt){
	matrix phi_1; phi_1.resize(N+2, 1); //phi at n+1th timestep
	double velocity;

	for (int i=1; i<N+1; i++){
		if (get_sgn(phi(i)) < 0){ //taking >= 0 as the other material
			velocity = var1->U(i+1, 1)/var1->U(i+1, 0); //"advection" velocity for the levelset function, phi
		}
		else {
			velocity = var2->U(i+1, 1)/var2->U(i+1, 0);
		}
		phi_1(i) = HJ_FirstOrder(velocity, dt, i);
	}

	phi = phi_1;
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::solver_MUSCL(EOS* eos1, EOS* eos2, gfmTests Test){
	
	MUSCL* m1 = dynamic_cast<MUSCL*>(var1);
	MUSCL* m2 = dynamic_cast<MUSCL*>(var2);
	slopeLimiter a;

	if(m1 && m2){
		a = m1->getLimiter();
	}
	else {
		std::cout << "dynamic casting to type MUSCL failed" << std::endl;
	}

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		for (int i=2; i<N+2; i++){
			if (get_sgn(phi(i-1)) < 0){
				ghost_boundary(var1, eos1, var2, eos2, i); //matrix 2 contains ghost points
			}

			else {
				ghost_boundary(var2, eos2, var1, eos1, i); //matrix 1 contains ghost points
			}
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//reconstruct data to piecewise linear representation
		m1->data_reconstruction(a);
		m2->data_reconstruction(a);

		//compute fluxes at current timestep		
		for (int i=1; i<N+2; i++){
			var1->compute_fluxes(eos1, i);
			var2->compute_fluxes(eos2, i);
		}

		//set timestep following CFL conditions with max wavespeed between both materials
		Smax = fmax(var1->Smax, var2->Smax);
		dt = CFL*(dx/Smax); 
		if (t + dt > Test.tstop) dt = Test.tstop - t;

		//evolving the levelset equation
		update_levelset_MUSCL(dt);

		//updating both materials, looping through real and ghost points
		for (int i=2; i<N+2; i++){
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

void GhostFluidMethods::solver_MUSCL(EOS* eos1, EOS* eos2, EOS* eos3, gfmTests Test){
	MUSCL* m1 = dynamic_cast<MUSCL*>(var1);
	MUSCL* m2 = dynamic_cast<MUSCL*>(var2);
	slopeLimiter a;

	if(m1 && m2){
		a = m1->getLimiter();
	}
	else {
		std::cout << "dynamic casting to type MUSCL failed" << std::endl;
	}

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){ //-1 or 0, the boundary is crossed
				if (get_sgn(phi(i)) < 0){//the interface is to the right of cell i
					for (int j=0; j<4; j++){
						//to the left of interface, the ghostfluid is in var2
						ghost_boundary(var1, eos1, var2, eos2, (i+1)-j); 
						//to the right of the interface, the ghostfluid is in var1
						ghost_boundary(var2, eos2, var1, eos1, (i+1)+j+1);
						//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
						//i, i-1 and i-2, while thst of var1 is i+1, i+2 and i+3.
					}
				}
				else {//the interface is to the left of cell i+1 (i+1 is negative while i is positive)
					for(int j=0; j<4; j++){
						//std::cout << var2->U.row(i+1-j) << std::endl;
						//to the left of interface, the ghostfluid is in var1
							//note that we are now in the right material which has EOS: eos3
						ghost_boundary(var2, eos2, var1, eos3, (i+1)-j);
						//to the right of the interface, the ghostfluid is in var2
						ghost_boundary(var1, eos3, var2, eos2, (i+1)+1+j);
						//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
						//i, i+1 and i+2, while thst of var1 is i-1, i-2 and i-3.
					}
				}
			} 
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//reconstruct data to piecewise linear representation
		m1->data_reconstruction(a);
		m2->data_reconstruction(a);

		//compute fluxes at current timestep		
		bool flag = true;
		//std::cout << "test1" << std::endl;
		for (int i=1; i<N+2; i++){
			if (phi(i-1) < 0 && flag==true){
				//std::cout <<"a " << i << std::endl;
				var1->compute_fluxes(eos1, i);
				//if (get_sgn(phi(i+1)) > 0) var1->compute_fluxes(eos1, i+1);
			}
			if (phi(i-1) >= 0){
				//std::cout <<"b " << i << std::endl;
				if (get_sgn(phi(i-2)) < 0) var2->compute_fluxes(eos1, i-1);
				var2->compute_fluxes(eos2, i);
				//if (get_sgn(phi(i+1)) < 0) var1->compute_fluxes(eos1, i+1);
				flag = false;
			}
			if (phi(i-1) < 0 && flag == false){
				//std::cout <<"c " << i << std::endl;
				if (get_sgn(phi(i-2)) > 0) var1->compute_fluxes(eos3, i-1);
				var1->compute_fluxes(eos3, i);
			}
			//if (count == 0 ) std::cout << var1->U.row(i) << std::endl; //'\t' << var2->U.row(i) << std::endl;
			//if (count ==1 )std::cout << i << '\t' << var1->F.row(i) << std::endl;//'\t' << var2->F.row(i) << std::endl;
		}

		//set timestep following CFL conditions with max wavespeed between both materials
		Smax = fmax(var1->Smax, var2->Smax);
		dt = CFL*(dx/Smax); 
		if (t + dt > Test.tstop) dt = Test.tstop - t;

		//updating both materials, looping through real and ghost points
		for (int i=2; i<N+2; i++){
			if (phi(i-1) < 0){
				var1->conservative_update_formula(dt, dx, i);
			}
			else {
				var2->conservative_update_formula(dt, dx, i);
			}
			//if (phi(i)<0) var1->conservative_update_formula(dt, dx, i);
			//else var2->conservative_update_formula(dt, dx, i);
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//evolving the levelset equation
		update_levelset_MUSCL(dt);
		
		t += dt;
		count += 1;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
}


void GhostFluidMethods::output_MUSCL(EOS* eos1, EOS* eos2){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	double d, u, P, e;
	for (int i=2; i<N+2; i++){
		if (get_sgn(phi(i-1)) < 0){
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
		outfile << X(i-1) << '\t' << d << '\t' << u
		<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done: GFM-MUSCL" << std::endl;
}

void GhostFluidMethods::output_MUSCL(EOS* eos1, EOS* eos2, EOS* eos3){
	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	double d, u, P, e;
	bool flag = true;
	for (int i=2; i<N+2; i++){
		if (phi(i-1) < 0 && flag == true){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos1->Pressure(var1->U, i);
			e = eos1->internalE(var1->U, i);
		}

		if (phi(i-1) >= 0) {
			d = var2->U(i, 0);
			u = var2->U(i, 1)/var2->U(i, 0);
			P = eos2->Pressure(var2->U, i);
			e = eos2->internalE(var2->U, i);
			flag = false;
		}
		if (phi(i-1) < 0 && flag == false){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos3->Pressure(var1->U, i);
			e = eos3->internalE(var1->U, i);
		}

		outfile << X(i) << '\t' << d << '\t' << u
		<< '\t' << P << '\t' << e << std::endl;
	}
}

//-----------------------------------------------------------------
//	EXACT
//-----------------------------------------------------------------

void GhostFluidMethods::exact_solver(gfmTests Test){
	//If there are 2 rarefraction waves, the exact solution no longer holds.
	double yL = Test.y1;
	double yR = Test.y2;

	vector WL = Test.initialL;
	vector WR = Test.initialR;

	//-------------------------------------------
	//constants
	double cL = sqrt(yL*WL(2)/WL(0));
	double cR = sqrt(yR*WR(2)/WR(0));

	double CONST1L = (yL-1)/(2*yL); // y-1 / 2y
	double CONST2L = (yL+1)/(2*yL); // y+1 / 2y
	double CONST3L = (2*yL)/(yL-1);
	double CONST4L = 2./(yL-1);
	double CONST5L = 2./(yL+1);
	double CONST6L = (yL-1)/(yL+1);
	double CONST7L = (yL-1)/2.;
	double CONST8L = yL-1;

	double CONST1R = (yR-1)/(2*yR); // y-1 / 2y
	double CONST2R = (yR+1)/(2*yR); // y+1 / 2y
	double CONST3R = (2*yR)/(yR-1);
	double CONST4R = 2./(yR-1);
	double CONST5R = 2./(yR+1);
	double CONST6R = (yR-1)/(yR+1);
	double CONST7R = (yR-1)/2.;

	double TOL = 1e-6;

	//data-dependent constants
	double AL = CONST5L/WL(0); 	double AR = CONST5R/WR(0);
	double BL = WL(2)*CONST6L; 	double BR = WR(2)*CONST6R;
	//-------------------------------------------

	auto check_pressure_pos_condition = [cL, cR, WL, WR, CONST4L, CONST4R]() { 
		//(Δu)crit ≡ 2aL/γ−1 + 2aR/γ−1 ≤ uR −uL , (4.82) Toro pg 127
		double du_crit = CONST4L*cL + CONST4R*cR;
		double du = WR(1) - WL(1);

		//If the pressure positivity condition is not satisfied, vacuum is generated
		// and the solver fails. Exit the program if this condition is not satisfied.
		if (du_crit <= du){ //This conditions ensures S*L <= S*R
			throw "Pressure positivity condition violated";
		}
	};

	//function of P, f(P, WL, WR)
	auto fL = [cL, AL, BL, WL, CONST4L, CONST1L](double P){
		double QL = sqrt(AL/(P + BL));
		double FL;

		if (P > WL(2)){ //Shock
			FL = (P - WL(2))*QL;
		}

		else { //Rarefraction
			FL = CONST4L*cL*(pow(P/WL(2), CONST1L) - 1);
		}

		return FL;
	};

	auto fR = [cR, AR, BR, WR, CONST4R, CONST1R](double P){
		double QR = sqrt(AR/(P + BR)); 
		double FR;

		if (P > WR(2)){ //Shock
			FR = (P - WR(2))*QR;
		}

		else { //Rarefraction
			FR = CONST4R*cR*(pow(P/WR(2), CONST1R) - 1);
		}

		return FR;
	};

	auto f = [WL, WR, fR, fL](double P){
		double du = (WR(1) - WL(1)); //velocity difference
		return fL(P) + fR(P) + du;   
	};

	auto fprimeL = [cL, AL, BL, WL, CONST2L](double P){
		double QL = sqrt(AL/(P + BL)); 
		double FLprime;
		if (P > WL(2)){ //Shock
			FLprime = QL*(1 - (P - WL(2))/(2.*(P + BL))); 
		}

		else { //Rarefraction
			FLprime = (1./(WL(0)*cL))*pow(P/WL(2), -CONST2L);
		}
		return FLprime;
	};

	auto fprimeR = [cR, AR, BR, WR, CONST2R](double P){
		double QR = sqrt(AR/(P + BR));
		double FRprime;
		if (P > WR(2)){ //Shock
			FRprime = QR*(1 - (P - WR(2))/(2.*(P + BR))); 
		}

		else { //Rarefraction
			FRprime = (1./(WR(0)*cR))*pow(P/WR(2), -CONST2R);
		}
		return FRprime;
	};

	auto fprime = [fprimeL, fprimeR](double P){
		return fprimeL(P) + fprimeR(P);  
	};

	auto newton_raphson = [f, fprime](double Pk){
		double Pk_1 = Pk - f(Pk)/fprime(Pk);
		return Pk_1;
	};

	auto relative_pressure_change = [](double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
		double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
		return CHA;
	};	

	auto compute_star_pressure = [yL, yR, cL, cR, AL, AR, BL, BR, WL, WR, CONST1L, CONST3L, CONST8L, 
		TOL, check_pressure_pos_condition, newton_raphson, relative_pressure_change](){
		//check positivity condition
		try {
			check_pressure_pos_condition();
		} 
		catch(const char* c){
			std::cout << c << std::endl;
			std::cout << "vacuum generated, terminating program" << std::endl;
		}

		double pPV, p0;
		
		//An approximation for p, p0 is required for the initial guess.
		//A poor choice of p0 results in the need for large number of iterations to achieve convergence
		
		pPV = 0.5*(WL(2) + WR(2)) - (1./8.)*(WR(1) - WL(1))*(WL(0) + WR(0))*(cL + cR);
		pPV = fmax(TOL, pPV);
		double Pmax = fmax(WL(2), WR(2));
		double Pmin = fmin(WL(2), WR(2));
		double Quser = 2.0;

		//if the pressure difference is small and the guess is within the initial pressure range
		if (Pmax/Pmin < Quser && pPV > Pmin && pPV < Pmax){
			//Select the linearised pressure guess
			p0 = pPV;
		}
		//If the guess is less than the initial pressure, two rarefraction wves have been formed
		else if (pPV < Pmin){
		//Since the exact solution for two rarefractions with a jump in y is invalid, calculation of
		//the pTR is only allowed for situations when yL = yR; 
			if (yL == yR){
				double pTR = pow((cL + cR - 0.5*CONST8L*(WR(1) - WL(1)))/((cL/pow(WL(2), CONST1L)) + (cR/pow(WR(2), CONST1L))), CONST3L);
				p0 = pTR;
			}
			else {
				throw "Exact solver attempted on two rarefraction waves with an EOS jump.";
			}
		}
		//If the pressure difference is large or if the guess is larger than Pmax
		else {
			//Select the two shock initial guess using pPV as estimate
			double QL = sqrt(AL/(pPV + BL));
			double QR = sqrt(AR/(pPV + BR));
			double pTS = (QL*WL(2) + QR*WR(2) - (WR(1) - WL(1)))/(QL + QR);
			p0 = pTS;
		}

		double Pk; Pk = p0; //First guess
		double Pk_1;
		double CHA = 0;
		double mixTOL;

		int count = 0;
		std::cout << p0 << std::endl;
		do{
			Pk_1 = newton_raphson(Pk);
			CHA = relative_pressure_change(Pk_1, Pk);
			Pk = Pk_1; //Set the iterate as the new guess
			//std::cout << CHA << '\t' << Pk << std::endl;
			count += 1;
			mixTOL = TOL*(1 + Pk);

			if (CHA < mixTOL) break;
			if (count == 20) std::cout << "Warning, maximum iterations reached for Newton's method" << std::endl;
		}while(count < 20);

		std::cout << count <<std::endl;
		return Pk;
	};

	auto compute_star_velocity = [WL, WR, fL, fR](double pstar){
		double ustar = 0.5*(WL(1) + WR(1)) + 0.5*(fR(pstar) - fL(pstar));
		return ustar;
	};

	auto compute_shock_density_L = [WL, CONST6L](double pstar){
		//turn this into a stored variable so the iteration does not have to be called multiple times
		//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
		double Pratio = pstar/WL(2);
		double dshock_L = WL(0)*((CONST6L + Pratio)/(CONST6L*Pratio + 1));
		return dshock_L;
	};

	auto compute_shock_density_R = [WR, CONST6R](double pstar){
		//turn this into a stored variable so the iteration does not have to be called multiple times
		//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
		double Pratio = pstar/WR(2);
		double dshock_R = WR(0)*((CONST6R + Pratio)/(CONST6R*Pratio + 1));
		return dshock_R;
	};

	auto compute_rarefraction_density_L = [yL, WL](double pstar){
		double drare_L = WL(0)*pow(pstar/WL(2), 1./yL);
		return drare_L;
	};

	auto compute_rarefraction_density_R = [yR, WR](double pstar){
		double drare_R= WR(0)*pow(pstar/WR(2), 1./yR);
		return drare_R;
	};

/*-----------------------------------------------------------
	Sampling
-----------------------------------------------------------*/
	//Solution of wavestructure within the domain of interest
	double pstar = compute_star_pressure();
	double ustar = compute_star_velocity(pstar);
	double dshockL = compute_shock_density_L(pstar);
	double dshockR = compute_shock_density_R(pstar);
	double drareL = compute_rarefraction_density_L(pstar);
	double drareR = compute_rarefraction_density_R(pstar);

	//primitive variables of left and right states
	double dL = WL(0); double dR = WR(0);
	double uL = WL(1); double uR = WR(1);
	double pL = WL(2); double pR = WR(2);

	//Shock Speeds
	double SL = uL - cL*sqrt((CONST2L*(pstar/pL) + CONST1L)); 
	double SR = uR + cR*sqrt((CONST2R*(pstar/pR) + CONST1R));

	//Rarefraction wave speeds
	double cstarL = cL*pow(pstar/pL, CONST1L);
	double cstarR = cR*pow(pstar/pR, CONST1R);

	double SHL = uL - cL; double SHR = uR + cR; //Head of fan
	double STL = ustar - cstarL; double STR = ustar + cstarR; //tail of fan
	/*
	std::cout << drareL << '\t' << drareR << std::endl;
	std::cout << SHL << '\t' << SHR << std::endl;
	std::cout << STL << '\t' << STR<< std::endl; */

	int newN = 1000;
	double newdx = Test.L/static_cast<double>(newN);
	matrix W(newN, 3);
	//Sampling is based on the position of wave with respect to time in x-t space,
	// characterised by its "speed" S = x/t.
	//When the solution at a specified time t is required the solution profiles are only a function of space x. Toro 137
	for (int i=0; i<newN; i++){
		double xPos = (static_cast<double>(i) + 0.5)*newdx;
		double S = (xPos - x0)/Test.tstop;
		//--------------------------------
		// Sampled region lies to the...
		//--------------------------------
		//Left side of Contact wave
		if (S <= ustar){
			//Left Shock wave
			if (pstar > pL){
				//Left of Shock Wave (unaffected by shock)
				if (S < SL){
					W.row(i) = WL;
				}
				//Right of Shock Wave (shocked material, star state)
				else {
					W(i, 0) = dshockL;
					W(i, 1) = ustar;
					W(i, 2) = pstar;
				}
			}
			//Left Rarefraction fan
			else {
				//Left of fastest rarefraction wave (unaffected by rarefraction)
				if (S < SHL){
					W.row(i) = WL;
				}
				//Out of the rarefraction fan, within the left star state
				else if (S > STL){
					W(i, 0) = drareL;
					W(i, 1) = ustar;
					W(i, 2) = pstar;
					//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
				}
				//Within the rarefraction fan
				else {
					double dLfan = dL*pow(CONST5L + (CONST6L/cL)*(uL - S), CONST4L);
					double uLfan = CONST5L*(cL + CONST7L*uL + S);
					double pLfan = pL*pow(CONST5L + (CONST6L/cL)*(uL - S), CONST3L);
					vector WLfan(dLfan, uLfan, pLfan);
					W.row(i) = WLfan;
				}
			}
		}
		//right side of contact wave
		else {
			//Right Shock wave
			if (pstar > pR){
				//Right of Shock Wave (unaffected by shock)
				if (S > SR){
					//std::cout << "Shock WR" << '\t' << i << std::endl;
					W.row(i) = WR;
				}
				//Left of Shock Wave (shocked material, star state)
				else {
					//std::cout << "Shock Star" << '\t' << i << std::endl;
					W(i, 0) = dshockR;
					W(i, 1) = ustar;
					W(i, 2) = pstar;
				}
			}
			//Right Rarefraction fan
			else {
				//Right of fastest rarefraction wave (unaffected by rarefraction)
				if (S > SHR){
					//std::cout << "rarefraction WR" << '\t' << i << std::endl;
					W.row(i) = WR;
				}
				//Out of the rarefraction fan, within the right star state
				else if (S < STR){
					//std::cout << "rarefraction Star" << '\t' << i << std::endl;
					W(i, 0) = drareR;
					W(i, 1) = ustar;
					W(i, 2) = pstar;
					//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
				}
				//Within the rarefraction fan
				else {
					//std::cout << "rarefraction fan" << '\t' << i << std::endl;
					double dRfan = dR*pow(CONST5R - (CONST6R/cR)*(uR - S), CONST4R);
					double uRfan = CONST5R*(-cR + CONST7R*uR + S);
					double pRfan = pR*pow(CONST5R - (CONST6R/cR)*(uR - S), CONST3R);
					vector WRfan(dRfan, uRfan, pRfan);
					W.row(i) = WRfan;
				}
			}
		}			
	}

	std::cout << compute_star_pressure() << std::endl;

	/*-------------------------------
		output
	-------------------------------*/
	std::ofstream outfile;
	outfile.open("dataexact.txt");

	for (int i=0; i<newN; i++){
		double xPos = (static_cast<double>(i) + 0.5)*newdx;
		double S = (xPos - x0)/Test.tstop;
		double e;
		if (S <= ustar){
			e = W(i, 2)/(W(i, 0)*(yL-1));
		}

		else {
			e = W(i, 2)/(W(i, 0)*(yR-1));
		}

		outfile << xPos << '\t' << W(i, 0) << '\t' << W(i, 1)
				<< '\t' << W(i, 2) << '\t' << e << std::endl;
	}
	outfile.close();
	std::cout << "done: exact" << std::endl;
}










