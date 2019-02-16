#include "GhostFluid.h"

GhostFluidMethods::GhostFluidMethods(double c, gfmTests Test)
	: LevelSetFunction(Test), CFL(c), Smax(0), dt(0), count(0), C(Test.number_of_materials, 11){
		var1 = new MUSCL(Test);
		var2 = new MUSCL(Test);
}

void GhostFluidMethods::ghost_boundary(MUSCL* Ureal, EOS* eosReal, MUSCL* Ughost, EOS* eosGhost, int i){

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

void GhostFluidMethods::ghost_boundary_RP(MUSCL* Uleft, EOS* eosleft, MUSCL* Uright, EOS* eosright, int i){
//Exact riemannsolver for idealgas and stiffened gas. Ghost point values are the RP star states.
//Consider an interface between i and i+1
	//converting to primitive variables
	double dl = Uleft->U(i-1, 0);
	double dr = Uright->U(i+2, 0);

	double ul = Uleft->U(i-1, 1)/Uleft->U(i-1, 0);
	double ur = Uright->U(i+2, 1)/Uright->U(i+2, 0);

	double Pl = eosleft->Pressure(Uleft->U, i-1);
	double Pr = eosright->Pressure(Uright->U, i+2);

	vector WL(dl, ul, Pl);
	vector WR(dr, ur, Pr);

	eosleft->y_constants(WL);
	eosright->y_constants(WR);

	//Solution of wavestructure within the domain of interest
	double pstar = compute_star_pressure(eosleft, eosright);
	//std::cout << "pstar = " << pstar << std::endl;
	double ustar = compute_star_velocity(pstar, eosleft, eosright);
	double dshockL = compute_shock_density(pstar, eosleft);
	double dshockR = compute_shock_density(pstar, eosright);
	double drareL = compute_rarefraction_density(pstar, eosleft);
	double drareR = compute_rarefraction_density(pstar, eosright);

	double Lghostdensity, Rghostdensity;
	if (pstar > Pl){
		//std::cout << "dshockL" << std::endl;
		Lghostdensity = dshockL;
	}
	else{
		//std::cout << "drareL" << std::endl;
		Lghostdensity = drareL;
	}

	if (pstar > Pr){
		//std::cout << "dshockR" << std::endl;
		Rghostdensity = dshockR;
	}
	else{
		//std::cout << "drareR" << std::endl;
		Rghostdensity = drareR;
	}
	//ghost fluid for left and right material
	vector WLghost(Lghostdensity, ustar, pstar);
	vector WRghost(Rghostdensity, ustar, pstar);

	//Redefining the real fluids at the boundary
	Uleft->U.row(i) = eosleft->conservedVar(WLghost);
	Uright->U.row(i+1) = eosright->conservedVar(WRghost);
	for (int j=0; j<4; j++){
		//Defining the states at ghost fluid cells
		Uleft->U.row(i+j+1) = eosleft->conservedVar(WLghost);
		Uright->U.row(i-j) = eosright->conservedVar(WRghost);
	}
}

//--------------------------------------------------------
//	MUSCL
//--------------------------------------------------------

void GhostFluidMethods::initial_conditions(EOS* eos1, EOS* eos2, gfmTests Test){
	vector euler1;
	vector euler2;

	//setting material EOS parameter
	eos1->y = Test.yL;
	eos2->y = Test.yR;

	//initial values for conserved variables
	euler1 = eos1->conservedVar(Test.initialL);
	euler2 = eos2->conservedVar(Test.initialR);


	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		signed_distance_function_1D(i); //refers to x(i+1)
		if (X(i+1)  < x0){ 
		//note, X(i+1) cant be set as <= x0, unless the level set is changed
		//to also allow this equality
			var1->U.row(i+2) = euler1; //setting real values of material 1
		}
		else{
			var2->U.row(i+2) = euler2; //setting real values of material 2
		}
	}

	for (int i=0; i<N+4; i++){
		//std::cout << var1->U.row(i) << std::endl;
	}

	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=1; i<N+1; i++){
		int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
		//std::cout << phi(i) << '\t' << testsgn << std::endl;
		if (testsgn > -2 && testsgn < 1 && i!=N){
			for (int j=0; j<3; j++){
				ghost_boundary(var1, eos1, var2, eos2, (i+1)-j); //matrix 2 contains ghost points
				ghost_boundary(var2, eos2, var1, eos1, (i+1)+j+1); //matrix 1 contains ghost points
			}
		}
	}


	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction

}

void GhostFluidMethods::initial_conditions(EOS* eos1, EOS* eos2, EOS* eos3, gfmTests Test){

	//setting material EOS parameter
	eos1->y = Test.yL;
	eos3->y = Test.yR;
	if (Test.yM1 != 0){
		eos2->y = Test.yM1;
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

	for (int i=1; i<N+1; i++){
		//std::cout << i << '\t' << get_sgn(phi(i)) << '\t' <<var1->U.row(i+1) << '\t' << var2->U.row(i+1) << std::endl;
	}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction

}

void GhostFluidMethods::initial_conditions(EOS* eos1, EOS* eos2, EOS* eos3, EOS* eos4, gfmTests Test){

	//setting material EOS parameter
	eos1->y = Test.yL;
	eos4->y = Test.yR;
	if (Test.yM1 != 0 && Test.yM2 != 0){
		eos2->y = Test.yM1;
		eos3->y = Test.yM2;
	}
	else {
		throw "Number of materials do not match the chosen solver.";
	}

	//initial values for conserved variables
	vector euler1 = eos1->conservedVar(Test.initialL);
	vector euler2 = eos2->conservedVar(Test.initialM1);
	vector euler3 = eos3->conservedVar(Test.initialM2);
	vector euler4 = eos4->conservedVar(Test.initialR);

	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		signed_distance_function_1D_3(i);
		if (X(i+1)  <= x0){
			var1->U.row(i+2) = euler1; //Left IC
		}
		else if (X(i+1) > x0 && X(i+1) <= x1){
			var2->U.row(i+2) = euler2; //Middle left IC
		}
		else if (X(i+1) > x1 && X(i+1) <= x2){
			var1->U.row(i+2) = euler3; //Middle right IC
		}
		else {
			var2->U.row(i+2) = euler4;
		}
		//std::cout << var1->U.row(i) << '\t' << var2->U.row(i) << std::endl;
		//std::cout << phi(i+1) << std::endl;
	}

	//////////////////////////////////////////////////////////
	//Populating the Ghost cells
	//////////////////////////////////////////////////////////
	bool flag1 = false; //Markers to identify which domain have been computed
	bool flag2 = false;
	bool flag3 = false;

	for (int i=1; i<N+1; i++){
		int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
		//std::cout << phi(i) << '\t' << testsgn << std::endl;
		if (testsgn > -2 && testsgn < 1 && i!=N){ //-1 or 0, the boundary is crossed
			if (get_sgn(phi(i)) < 0){//the interface is to the right of cell i
				if (flag1 == false && flag2 == false && flag3 == false){
					for (int j=0; j<4; j++){
						//std::cout << "flag1" <<std::endl;
						//to the left of interface, the ghostfluid is in var2
						ghost_boundary(var1, eos1, var2, eos2, (i+1)-j); 
						//to the right of the interface, the ghostfluid is in var1
						ghost_boundary(var2, eos2, var1, eos1, (i+1)+j+1);
						//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
						//i, i-1 and i-2, while thst of var1 is i+1, i+2 and i+3.
						if (j >= 3) flag1 = true;
					}
				}
				//populating the 3rd boundary, second discontinuity in var1
				//the materials now obey eos3 in var1 and eos4 in var2
				else if (flag1 == true && flag2 == true && flag3 == false){
					for (int j=0; j<4; j++){
						//std::cout << "flag3" <<std::endl;
						ghost_boundary(var1, eos3, var2, eos4, (i+1)-j); 
						ghost_boundary(var2, eos4, var1, eos3, (i+1)+j+1);
						if (j >= 3) flag3 = true;
					}
				}
			}
			else {//the interface is to the left of cell i+1 (i+1 is negative while i is positive)
				for(int j=0; j<4; j++){
					//to the left of interface, the ghostfluid is in var1
					//std::cout << var2->U.row(i+1-j) << std::endl;
					if(flag1 == true && flag2 == false && flag3 == false){
						//std::cout << "flag2" <<std::endl;
						//note that we are now in the right material which has EOS: eos3
						ghost_boundary(var2, eos2, var1, eos3, (i+1)-j);
						//to the right of the interface, the ghostfluid is in var2
						ghost_boundary(var1, eos3, var2, eos2, (i+1)+1+j);
						//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
						//i, i+1 and i+2, while thst of var1 is i-1, i-2 and i-3.
						if (j >= 3) flag2 = true;
					}
				}
			}
		}
	}

	//for (int i=1; i<N+1; i++){
	//	std::cout << i << '\t' << get_sgn(phi(i)) << '\t' <<var1->U.row(i+1) << '\t' << var2->U.row(i+1) << std::endl;
	//}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::update_levelset(double dt){
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

void GhostFluidMethods::solver(EOS* eos1, EOS* eos2, gfmTests Test){
	
	slopeLimiter a;

	a = var1->getLimiter();

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){
				for (int j=0; j<3; j++){
				ghost_boundary(var1, eos1, var2, eos2, (i+1)-j); //matrix 2 contains ghost points
				ghost_boundary(var2, eos2, var1, eos1, (i+1)+j+1); //matrix 1 contains ghost points
				}
			}
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//reconstruct data to piecewise linear representation
		var1->data_reconstruction(a);
		var2->data_reconstruction(a);

		//compute fluxes at current timestep		
		for (int i=1; i<N+2; i++){
			//if (count==0)std::cout << i << '\t' << phi(i-1) << std::endl;
			if(phi(i-1) < 0){
				//std::cout << "phi < 0" << std::endl;
				var1->compute_fluxes(eos1, i);
			}
			if (phi(i-1) > 0){
				//std::cout << "phi > 0" << std::endl;
				//std::cout << i << std::endl;
				if (get_sgn(phi(i-2)) <= 0) {
					//std::cout << "phi > 0 and phi-1 < 0" << std::endl;
					var2->compute_fluxes(eos2, i-1);
				}
				var2->compute_fluxes(eos2, i);
			}
			//if (count==0)std::cout << i << '\t' << var1->U.row(i+1) << '\t' << '\t' << var2->U.row(i+1) << std::endl;
			//if (count==0)std::cout << i << '\t' << var1->F.row(i) << '\t' << '\t' << var2->F.row(i) << std::endl;
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
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//evolving the levelset equation
		update_levelset(dt);
		
		t += dt;
		count += 1;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
}

void GhostFluidMethods::solver(EOS* eos1, EOS* eos2, EOS* eos3, gfmTests Test){

	slopeLimiter a;

	a = var1->getLimiter();

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
		var1->data_reconstruction(a);
		var2->data_reconstruction(a);

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
		update_levelset(dt);
		
		t += dt;
		count += 1;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
}


void GhostFluidMethods::solver(EOS* eos1, EOS* eos2, EOS* eos3, EOS* eos4, gfmTests Test){
	slopeLimiter a;

	a = var1->getLimiter();

	double t = 0.0;
	do{
		//compute ghost fluid boundaries
		bool flag1 = false; //Markers to identify which domain have been computed
		bool flag2 = false;
		bool flag3 = false;

		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){ //-1 or 0, the boundary is crossed
				if (get_sgn(phi(i)) < 0){//the interface is to the right of cell i
					if (flag1 == false && flag2 == false && flag3 == false){
						for (int j=0; j<4; j++){
							//std::cout << "flag1" <<std::endl;
							//to the left of interface, the ghostfluid is in var2
							ghost_boundary(var1, eos1, var2, eos2, (i+1)-j); 
							//to the right of the interface, the ghostfluid is in var1
							ghost_boundary(var2, eos2, var1, eos1, (i+1)+j+1);
							//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
							//i, i-1 and i-2, while thst of var1 is i+1, i+2 and i+3.
							if (j >= 3) flag1 = true;
						}
					}
					//populating the 3rd boundary, second discontinuity in var1
					//the materials now obey eos3 in var1 and eos4 in var2
					else if (flag1 == true && flag2 == true && flag3 == false){
						for (int j=0; j<4; j++){
							//std::cout << "flag3" <<std::endl;
							ghost_boundary(var1, eos3, var2, eos4, (i+1)-j); 
							ghost_boundary(var2, eos4, var1, eos3, (i+1)+j+1);
							if (j >= 3) flag3 = true;
						}
					}
				}
				else {//the interface is to the left of cell i+1 (i+1 is negative while i is positive)
					for(int j=0; j<4; j++){
						//to the left of interface, the ghostfluid is in var1
						//std::cout << var2->U.row(i+1-j) << std::endl;
						if(flag1 == true && flag2 == false && flag3 == false){
							//std::cout << "flag2" <<std::endl;
							//note that we are now in the right material which has EOS: eos3
							ghost_boundary(var2, eos2, var1, eos3, (i+1)-j);
							//to the right of the interface, the ghostfluid is in var2
							ghost_boundary(var1, eos3, var2, eos2, (i+1)+1+j);
							//since cell i lies in the real fluid of var 1, ghost fluid of var2 goes from
							//i, i+1 and i+2, while thst of var1 is i-1, i-2 and i-3.
							if (j >= 3) flag2 = true;
						}
					}
				}
			}
		}

		if(flag1==false || flag2==false || flag3==false){
			throw "Failed to populate Ghost Fluid Cells";
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		//reconstruct data to piecewise linear representation
		var1->data_reconstruction(a);
		var2->data_reconstruction(a);

		//compute fluxes at current timestep		
		//using the previous flags which should now be true,
		for (int i=1; i<N+2; i++){
			if (phi(i-1) < 0 && flag1==true && flag2==true){
				//std::cout <<"a " << i << std::endl;
				var1->compute_fluxes(eos1, i);
				//if (get_sgn(phi(i+1)) > 0) var1->compute_fluxes(eos1, i+1);
			}
			if (phi(i-1) >= 0 && flag2==true){
				//std::cout <<"b " << i << std::endl;
				if (get_sgn(phi(i-2)) < 0) var2->compute_fluxes(eos1, i-1);
				var2->compute_fluxes(eos2, i);
				//if (get_sgn(phi(i+1)) < 0) var1->compute_fluxes(eos1, i+1);
				flag1 = false;
			}
			if (phi(i-1) < 0 && flag1 == false){
				//std::cout <<"c " << i << std::endl;
				if (get_sgn(phi(i-2)) > 0) var1->compute_fluxes(eos3, i-1);
				var1->compute_fluxes(eos3, i);
				flag2 = false;
			}
			if (phi(i-1) >= 0 && flag1==false && flag2==false){
				if (get_sgn(phi(i-2)) < 0) var2->compute_fluxes(eos1, i-1);
				var2->compute_fluxes(eos2, i);
			}
			//if (count == 1 ) std::cout << var1->U.row(i) << std::endl; //'\t' << var2->U.row(i) << std::endl;
			//if (count == 180 )std::cout << i << '\t' << var1->F.row(i) << std::endl;//'\t' << var2->F.row(i) << std::endl;
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
		update_levelset(dt);
		
		t += dt;
		count += 1;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
}

///////////////////////////////////////////////////////////////////
// Riemann Problem based GFM
///////////////////////////////////////////////////////////////////

void GhostFluidMethods::initial_conditions_RP(EOS* eos1, EOS* eos2, gfmTests Test){

	vector euler1;
	vector euler2;

	//setting material EOS parameter
	eos1->y = Test.yL;
	eos2->y = Test.yR;

	eos1->y_constants(Test.initialL);
	eos2->y_constants(Test.initialR);

	//initial values for conserved variables
	euler1 = eos1->conservedVar(Test.initialL);
	euler2 = eos2->conservedVar(Test.initialR);


	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		signed_distance_function_1D(i); //refers to x(i+1)
		if (X(i+1)  < x0){
			var1->U.row(i+2) = euler1; //setting real values of material 1
		}
		else{
			var2->U.row(i+2) = euler2; //setting real values of material 2
		}
	}


	//After setting the real fluids, populate the empty cells with ghost values
	for (int i=1; i<N+1; i++){
		int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
		//std::cout << phi(i) << '\t' << testsgn << std::endl;
		if (testsgn > -2 && testsgn < 1 && i!=N){
			//std::cout << "boundary is at i = " << i << std::endl;
			ghost_boundary_RP(var1, eos1, var2, eos2, i+1); //matrix 2 contains ghost points
		}
	}

	for (int i=0; i<N; i++){
		//std::cout << i << '\t' << phi(i+1) << '\t' << get_sgn(phi(i+1)) << var1->U.row(i+2) << '\t' << '\t' << var2->U.row(i+2) << std::endl;
	}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::solver_RP(EOS* eos1, EOS* eos2, gfmTests Test){
	slopeLimiter a;

	a = var1->getLimiter();

	double t = 0.0;

	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){
				//std::cout << count << " boundary is at i = " << i << std::endl;
				ghost_boundary_RP(var1, eos1, var2, eos2, i+1); //matrix 2 contains ghost points
			}

			//if (count==0) std::cout << var1->U.row(i) << std::endl; 
		}

		var1->boundary_conditions();
		var2->boundary_conditions();

		for (int i=2; i<N+2; i++){
			//std::cout << var1->U.row(i) << std::endl;
		}

		//reconstruct data to piecewise linear representation
		var1->data_reconstruction(a);
		var2->data_reconstruction(a);


		for (int i=1; i<N+2; i++){
			//if (count==0)std::cout << i << '\t' << phi(i-1) << std::endl;
			if(phi(i-1) < 0){
				var1->compute_fluxes(eos1, i);
			}
			if (phi(i-1) >= 0){
				var2->compute_fluxes(eos2, i);
			}
			//if (count==0)std::cout << i << '\t' << var1->U.row(i+1) << '\t' << '\t' << var2->U.row(i+1) << std::endl;
			//if (count==1)std::cout << i << '\t' << var1->F.row(i) << '\t' << '\t' << var2->F.row(i) << std::endl;
		}

		for (int i=1; i<N+2; i++){
			//if (count == 1)std::cout << i << '\t' << var1->F.row(i) << '\t' << '\t' << '\t' << var2->F.row(i) << std::endl;
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
			if (phi(i-1) >=0) {
				var2->conservative_update_formula(dt, dx, i);
			}
			//if (count==0) std::cout << var1->U.row(i) << std::endl;
		}

		//evolving the levelset equation
		update_levelset(dt);

		var1->boundary_conditions();
		var2->boundary_conditions();
		
		t += dt;
		count += 1;

	}while(t<Test.tstop);
	std::cout << count << std::endl;

}
///////////////////////////////////////////////////////////////////


void GhostFluidMethods::output(EOS* eos1, EOS* eos2){

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

void GhostFluidMethods::output(EOS* eos1, EOS* eos2, EOS* eos3){
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

		outfile << X(i-1) << '\t' << d << '\t' << u
		<< '\t' << P << '\t' << e << std::endl;
	}
}

void GhostFluidMethods::output(EOS* eos1, EOS* eos2, EOS* eos3, EOS* eos4){
	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	double d, u, P, e;
	bool flag1 = true;
	bool flag2 = true;

	for (int i=2; i<N+2; i++){
		if (phi(i-1) < 0 && flag1==true && flag2==true){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos1->Pressure(var1->U, i);
			//e = eos1->internalE(var1->U, i);
			e = eos1->internalE(var1->U, i)*var1->U(i, 0);
		}

		if (phi(i-1) >= 0 && flag2==true) {
			d = var2->U(i, 0);
			u = var2->U(i, 1)/var2->U(i, 0);
			P = eos2->Pressure(var2->U, i);
			//e = eos2->internalE(var2->U, i);
			e = eos2->internalE(var2->U, i)*var2->U(i, 0);
			flag1 = false;
		}
		if (phi(i-1) < 0 && flag1 == false){
			d = var1->U(i, 0);
			u = var1->U(i, 1)/var1->U(i, 0);
			P = eos3->Pressure(var1->U, i);
			//e = eos3->internalE(var1->U, i);
			e = eos3->internalE(var1->U, i)*var1->U(i, 0);
			flag2 = false;
		}
		if (phi(i-1) >= 0 && flag1==false && flag2==false){
			d = var2->U(i, 0);
			u = var2->U(i, 1)/var2->U(i, 0);
			P = eos4->Pressure(var2->U, i);
			//e = eos4->internalE(var2->U, i);
			e = eos4->internalE(var2->U, i)*var2->U(i, 0);
		}

		outfile << X(i-1) << '\t' << d << '\t' << u
		<< '\t' << P << '\t' << e << std::endl;
	}
}

///////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------
//	EXACT
//-----------------------------------------------------------------
///////////////////////////////////////////////////////////////////

double GhostFluidMethods::f(double P, EOS* eosleft, EOS* eosright){
	double du = (eosright->C(12) - eosleft->C(12)); //velocity difference
	return eosleft->fk(P) + eosright->fk(P) + du;   
}

double GhostFluidMethods::fprime(double P, EOS* eosleft, EOS* eosright){
	return eosleft->fprimek(P) + eosright->fprimek(P);   
}

void GhostFluidMethods::check_pressure_pos_condition(EOS* eosleft, EOS* eosright) { 
	//(Δu)crit ≡ 2aL/γ−1 + 2aR/γ−1 ≤ uR −uL , (4.82) Toro pg 127
	double du_crit = eosleft->C(4)*eosleft->C(0) + eosright->C(4)*eosright->C(0);
	double du = eosright->C(12) - eosleft->C(12);

	//If the pressure positivity condition is not satisfied, vacuum is generated
	// and the solver fails. Exit the program if this condition is not satisfied.
	if (du_crit <= du){ //This conditions ensures S*L <= S*R
		throw "Pressure positivity condition violated";
	}
}

double GhostFluidMethods::newton_raphson(double Pk, EOS* eosleft, EOS* eosright){
	double Pk_1 = Pk - f(Pk, eosleft, eosright)/fprime(Pk, eosleft, eosright);
	return Pk_1;
}

double GhostFluidMethods::compute_star_pressure(EOS* eosleft, EOS* eosright){
//------------------------------------------
	double TOL = 1e-6;
	auto relative_pressure_change = [](double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
		double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
		return CHA;
	};	
//-------------------------------------------
	//check positivity condition
	try {
		check_pressure_pos_condition(eosleft, eosright);
	} 
	catch(const char* c){
		std::cout << c << std::endl;
		std::cout << "vacuum generated, terminating program" << std::endl;
	}
	double pPV, p0;
	
	//An approximation for p, p0 is required for the initial guess.
	//A poor choice of p0 results in the need for large number of iterations to achieve convergence
	double soundspeed1 = eosleft->C(0);
	double soundspeed2 = eosright->C(0);
	pPV = 0.5*(eosleft->C(13) + eosright->C(13)) - (1./8.)*(eosright->C(12) - eosleft->C(12))*(eosleft->C(11) + eosright->C(11))*(soundspeed1 + soundspeed2);
	pPV = fmax(TOL, pPV);
	double Pmax = fmax(eosleft->C(13), eosright->C(13));
	double Pmin = fmin(eosleft->C(13), eosright->C(13));
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
		double y1 = eosleft->C(10); double y2 = eosright->C(10);
		if (y1 == y2){
			double pTR = pow((soundspeed1 + soundspeed2 - 0.5*(y1-1)*(eosright->C(12) - eosleft->C(12)))/((soundspeed1/pow(eosleft->C(13), eosleft->C(1))) + (soundspeed2/pow(eosright->C(13), eosright->C(1)))), eosleft->C(3));
			p0 = pTR;
		}
		else {
			throw "Exact solver attempted on two rarefraction waves with an EOS jump.";
		}
	}
	//If the pressure difference is large or if the guess is larger than Pmax
	else {
		//Select the two shock initial guess using pPV as estimate
		double Q1 = sqrt(eosleft->C(8)/(pPV + eosleft->C(9)));
		double Q2 = sqrt(eosright->C(8)/(pPV + eosright->C(9)));
		double pTS = (Q1*eosleft->C(13) + Q2*eosright->C(13) - (eosright->C(12) - eosright->C(12)))/(Q1 + Q2);
		p0 = pTS;
	}

	std::cout << "initial pressure = "<< p0 << std::endl;
	double Pk; Pk = p0; //First guess
	double Pk_1;
	double CHA = 0;

	int count = 0;
	//std::cout << p0 << std::endl;
	do{
		Pk_1 = newton_raphson(Pk, eosleft, eosright);
		CHA = relative_pressure_change(Pk_1, Pk);
		Pk = Pk_1; //Set the iterate as the new guess
		//std::cout << CHA << '\t' << Pk << std::endl;
		count += 1;

		if (CHA < TOL) break;
		if (count == 20) std::cout << "Warning, maximum iterations reached for Newton's method" << std::endl;
	}while(count < 20);

	//std::cout << count <<std::endl;
	return Pk;
}

double GhostFluidMethods::compute_star_velocity(double pstar, EOS* eosleft, EOS* eosright){
	double ustar = 0.5*(eosleft->C(12) + eosright->C(12)) + 0.5*(eosright->fk(pstar) - eosleft->fk(pstar));
	return ustar;
}

double GhostFluidMethods::compute_shock_density(double pstar, EOS* eos){
	//turn this into a stored variable so the iteration does not have to be called multiple times
	//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
	double Pratio = pstar/eos->C(13);
	double dshock = eos->C(11)*((eos->C(6) + Pratio)/(eos->C(6)*Pratio + 1));
	return dshock;
}

double GhostFluidMethods::compute_rarefraction_density(double pstar, EOS* eos){
	double drare = eos->C(11)*pow(pstar/eos->C(13), 1./eos->C(10));
	return drare;
}

void GhostFluidMethods::exact_solver(gfmTests Test, EOS* eosleft, EOS* eosright){
	//If there are 2 rarefraction waves, the exact solution no longer holds.

	//setting material EOS parameter
	//eosleft->y = Test.yL;
	//eosright->y = Test.yR;

	vector WL = Test.initialL;	
	vector WR = Test.initialR;

	eosleft->y_constants(WL);
	eosright->y_constants(WR);

	double yL = eosleft->C(10);
	double yR = eosright->C(10);

	//-----------------------------------------------------------
	//Let the exact solution have a resolution of 1000 grid points
	int newN = 1000;
	double newdx = Test.L/static_cast<double>(newN);
	matrix W(newN, 3);
	//-----------------------------------------------------------

	if (Test.number_of_materials == 2){

		//-----------------Sampling------------------
		//Solution of wavestructure within the domain of interest
		double pstar = compute_star_pressure(eosleft, eosright);
		double ustar = compute_star_velocity(pstar, eosleft, eosright);
		double dshockL = compute_shock_density(pstar, eosleft);
		double dshockR = compute_shock_density(pstar, eosright);
		double drareL = compute_rarefraction_density(pstar, eosleft);
		double drareR = compute_rarefraction_density(pstar, eosright);

		//primitive variables of left and right states
		double dL = WL(0); double dR = WR(0);
		double uL = WL(1); double uR = WR(1);
		double pL = WL(2); double pR = WR(2);

		//Shock Speeds
		double SL = uL - eosleft->C(0)*sqrt((eosleft->C(2)*(pstar/pL) + eosleft->C(1))); 
		double SR = uR + eosright->C(0)*sqrt((eosright->C(2)*(pstar/pR) + eosright->C(1)));

		//Rarefraction wave speeds
		double cstarL = eosleft->C(0)*pow(pstar/pL, eosleft->C(1));
		double cstarR = eosright->C(0)*pow(pstar/pR, eosright->C(1));

		double SHL = uL - eosleft->C(0); double SHR = uR + eosright->C(0); //Head of fan
		double STL = ustar - cstarL; double STR = ustar + cstarR; //tail of fan
		
		//std::cout << drareL << '\t' << drareR << std::endl;
		//std::cout << "Shockspeeds = ";
		//std::cout << SL << '\t' << SR << std::endl;
		//std::cout << SHL << '\t' << SHR << std::endl;
		//std::cout << STL << '\t' << STR << std::endl; 

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
						double dLfan = dL*pow(eosleft->C(5) + (eosleft->C(6)/eosleft->C(0))*(uL - S), eosleft->C(4));
						double uLfan = eosleft->C(5)*(eosleft->C(0) + eosleft->C(7)*uL + S);
						double pLfan = pL*pow(eosleft->C(5) + (eosleft->C(6)/eosleft->C(0))*(uL - S), eosleft->C(3));
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
						double dRfan = dR*pow(eosright->C(5) - (eosright->C(6)/eosright->C(0))*(uR - S), eosright->C(4));
						double uRfan = eosright->C(5)*(-eosright->C(0) + eosright->C(7)*uR + S);
						double pRfan = pR*pow(eosright->C(5) - (eosright->C(6)/eosright->C(0))*(uR - S), eosright->C(3));
						vector WRfan(dRfan, uRfan, pRfan);
						W.row(i) = WRfan;
					}
				}
			}			
		}

		std::cout << "Pstar = " << pstar << '\t' << "ustar = " << ustar << std::endl;

		//-------------------------------
		//	output
		//-------------------------------
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

	else{
		throw "Exact solver is only avaliable for a single material interface.";
	}

}

void GhostFluidMethods::exact_solver(gfmTests Test, EOS* eosleft, EOS* eosmiddle, EOS* eosright){
//NOTE: only works on Fedkiw (2002), Test B and C
	//For 3 masterisls, Consider 2 RP between Left and Middle states, as well as Right and Middle states.
	//Assuming that the stop time is before interactions between the generated waves

	//setting material EOS parameter
	eosleft->y = Test.yL;
	eosright->y = Test.yR; //get rid of this!!!
	eosmiddle->y = Test.yM1;

	vector WL = Test.initialL;	
	vector WR = Test.initialR;
	vector WM1 = Test.initialM1;

	eosleft->y_constants(WL);
	eosright->y_constants(WR);
	eosmiddle->y_constants(WM1);
	
	double yL = eosleft->C(10);
	double yR = eosright->C(10);
	double yM1 = eosmiddle->C(10);

	//----------------Rename variables--------------------------
	//variablename_1 represents solution to the leftmost RP
	//----------------------------------------------------------
	//Solution of wavestructure within the domain of interest between WL and WM1
	//This is shared with the solver for 4 initial conditions
	double pstar_1 = compute_star_pressure(eosleft, eosmiddle);
	double ustar_1 = compute_star_velocity(pstar_1, eosleft, eosmiddle);
	double dshockL_1 = compute_shock_density(pstar_1, eosleft);
	double dshockR_1 = compute_shock_density(pstar_1, eosmiddle);
	double drareL_1 = compute_rarefraction_density(pstar_1, eosleft);
	double drareR_1 = compute_rarefraction_density(pstar_1, eosmiddle);

	//primitive variables of left and right states
	double dL_1 = WL(0); double dR_1 = WM1(0);
	double uL_1 = WL(1); double uR_1 = WM1(1);
	double pL_1 = WL(2); double pR_1 = WM1(2);

	//Shock Speeds
	double SL_1 = uL_1 - eosleft->C(0)*sqrt((eosleft->C(2)*(pstar_1/pL_1) + eosleft->C(1))); 
	double SR_1 = uR_1 + eosmiddle->C(0)*sqrt((eosmiddle->C(2)*(pstar_1/pR_1) + eosmiddle->C(1)));

	//Rarefraction wave speeds
	double cstarL_1 = eosleft->C(0)*pow(pstar_1/pL_1, eosleft->C(1));
	double cstarR_1 = eosmiddle->C(0)*pow(pstar_1/pR_1, eosmiddle->C(1));

	double SHL_1 = uL_1 - eosleft->C(0); double SHR_1 = uR_1 + eosmiddle->C(0); //Head of fan
	double STL_1 = ustar_1 - cstarL_1; double STR_1 = ustar_1 + cstarR_1; //tail of fan

//Tracking the right moving shock wave. For this problem, we know that the right moving wave of the
	//first riemann problem is a shock wave.
	//time when shock wave reaches the discontinuity
	double t_shock = Test.tstop - (Test.x1 - Test.x0)/SR_1;
	vector WstarL(dshockR_1, ustar_1, pstar_1);	
	//reassign eosleft with star state
	eosleft->y_constants(WstarL);

	//Solution of wavestructure within the domain of interest
	double pstar = compute_star_pressure(eosleft, eosright);
	double ustar = compute_star_velocity(pstar, eosleft, eosright);
	double dshockL = compute_shock_density(pstar, eosleft);
	double dshockR = compute_shock_density(pstar, eosright);
	double drareL = compute_rarefraction_density(pstar, eosleft);
	double drareR = compute_rarefraction_density(pstar, eosright);

	//primitive variables of left and right states
	double dL = WstarL(0); double dR = WR(0);
	double uL = WstarL(1); double uR = WR(1);
	double pL = WstarL(2); double pR = WR(2);

	double soundspeedL = eosleft->C(0);
	double soundspeedR = eosright->C(0);

	//Shock Speeds
	double SL = uL - soundspeedL*sqrt((eosmiddle->C(2)*(pstar/pL) + eosmiddle->C(1))); 
	double SR = uR + soundspeedR*sqrt((eosright->C(2)*(pstar/pR) + eosright->C(1)));

	//Rarefraction wave speeds
	double cstarL = soundspeedL*pow(pstar/pL, eosmiddle->C(1));
	double cstarR = soundspeedR*pow(pstar/pR, eosright->C(1));

	double SHL = uL - soundspeedL; double SHR = uR + soundspeedR; //Head of fan
	double STL = ustar - cstarL; double STR = ustar + cstarR; //tail of fan


	//-----------------------------------------------------------
	//	Samp[ling]
	//-----------------------------------------------------------
	//-----------------------------------------------------------
	//Let the exact solution have a resolution of 1000 grid points
	int newN = 1000;
	double newdx = Test.L/static_cast<double>(newN);
	matrix W(newN, 3);
	//-----------------------------------------------------------
	//consider the time after the wave has hit the discontinuity, forming a new riemann problem
	for (int i=0; i<newN; i++){

		double xPos = (static_cast<double>(i) + 0.5)*newdx;

		if (Test.tstop >= t_shock){
			double S = (xPos - Test.x1)/t_shock;
			if (S <= ustar){ //left of new contact wave
				//Left Shock wave
				if (pstar > pL){
					//Left of Shock Wave (unaffected by new shock)
					if (S < SL){
						W.row(i) = WstarL;
						//std::cout << "Shock WstarL" << '\t' << i << std::endl;
					}
					//Right of Shock Wave (shocked material, star state)
					else {
						W(i, 0) = dshockL;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
						//std::cout << "Shock new WstarL" << '\t' << i << std::endl;
					}
				}
				//Left Rarefraction fan
				else {
					//Left of fastest rarefraction wave (unaffected by rarefraction)
					if (S < SHL){
						W.row(i) = WstarL;
						//std::cout << "Rarefraction WstarL" << '\t' << i << std::endl;
					}
					//Out of the rarefraction fan, within the left star state
					else if (S > STL){
						W(i, 0) = drareL;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
						//std::cout << "Rarefraction new WstarL" << '\t' << i << std::endl;
					}
					//Within the rarefraction fan
					else {
						double dLfan = dL*pow(eosmiddle->C(5) + (eosmiddle->C(6)/soundspeedL)*(uL - S), eosmiddle->C(4));
						double uLfan = eosmiddle->C(5)*(soundspeedL + eosmiddle->C(7)*uL + S);
						double pLfan = pL*pow(eosmiddle->C(5) + (eosmiddle->C(6)/soundspeedL)*(uL - S), eosmiddle->C(3));
						vector WLfan(dLfan, uLfan, pLfan);
						W.row(i) = WLfan;
						//std::cout << "Rarefraction new WfanL" << '\t' << i << std::endl;
					}
				}
			}

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
						//std::cout << "Shock new WR" << '\t' << i << std::endl;
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
						//std::cout << "rarefraction new WR" << '\t' << i << std::endl;
						W(i, 0) = drareR;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
						//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
					}
					//Within the rarefraction fan
					else {
						//std::cout << "rarefraction WR" << '\t' << i << std::endl;
						double dRfan = dR*pow(eosright->C(5) - (eosright->C(6)/soundspeedR)*(uR - S), eosright->C(4));
						double uRfan = eosright->C(5)*(-soundspeedR + eosright->C(7)*uR + S);
						double pRfan = pR*pow(eosright->C(5) - (eosright->C(6)/soundspeedR)*(uR - S), eosright->C(3));
						vector WRfan(dRfan, uRfan, pRfan);
						W.row(i) = WRfan;
					}
				}
			}
		}

//----------//If the shockwave has not reached the discontinuity------------------
		else {
			double S = (xPos - Test.x0)/Test.tstop;
			if (S <= ustar_1){
				//Left Shock wave
				if (pstar_1 > pL_1){
					//Left of Shock Wave (unaffected by shock)
					if (S < SL_1){
						W.row(i) = WL;
						//std::cout << "Shock WL" << '\t' << i << std::endl;
					}
					//Right of Shock Wave (shocked material, star state)
					else {
						W(i, 0) = dshockL_1;
						W(i, 1) = ustar_1;
						W(i, 2) = pstar_1;
						//std::cout << "Shock WL star" << '\t' << i << std::endl;
					}
				}
				//Left Rarefraction fan
				else {
					//Left of fastest rarefraction wave (unaffected by rarefraction)
					if (S < SHL_1){
						W.row(i) = WL;
						//std::cout << "Rarefraction WL" << '\t' << i << std::endl;
					}
					//Out of the rarefraction fan, within the left star state
					else if (S > STL_1){
						W(i, 0) = drareL_1;
						W(i, 1) = ustar_1;
						W(i, 2) = pstar_1;
						//std::cout << "Rarefraction WL star" << '\t' << i << std::endl;
					}
					//Within the rarefraction fan
					else {
						double dLfan = dL_1*pow(eosleft->C(5) + (eosleft->C(6)/eosleft->C(0))*(uL_1 - S), eosleft->C(4));
						double uLfan = eosleft->C(5)*(eosleft->C(0) + eosleft->C(7)*uL_1 + S);
						double pLfan = pL_1*pow(eosleft->C(5) + (eosleft->C(6)/eosleft->C(0))*(uL_1 - S), eosleft->C(3));
						vector WLfan(dLfan, uLfan, pLfan);
						W.row(i) = WLfan;
						//std::cout << "Rarefraction WL fan" << '\t' << i << std::endl;
					}
				}
			}

			else {
				//Right Shock wave
				if (pstar_1 > pR_1){
					//Right of Shock Wave (unaffected by shock)
					if (S > SR_1){
						//std::cout << "Shock WM" << '\t' << i << std::endl;
						W.row(i) = WM1;

					}
					//Left of Shock Wave (shocked material, star state)
					else {
						//std::cout << "Shock WM Star" << '\t' << i << std::endl;
						W(i, 0) = dshockR_1;
						W(i, 1) = ustar_1;
						W(i, 2) = pstar_1;

					}
				}
				//Right Rarefraction fan
				else {
					//Right of fastest rarefraction wave (unaffected by rarefraction)
					if (S > SHR_1){
						//std::cout << "rarefraction WM" << '\t' << i << std::endl;
						W.row(i) = WM1;
					}
					//Out of the rarefraction fan, within the right star state
					else if (S < STR_1){
						//std::cout << "rarefraction WM Star" << '\t' << i << std::endl;
						W(i, 0) = drareR_1;
						W(i, 1) = ustar_1;
						W(i, 2) = pstar_1;
						//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
					}
					//Within the rarefraction fan
					else {
						//std::cout << "rarefraction WM fan" << '\t' << i << std::endl;
						double dRfan = dR_1*pow(eosmiddle->C(5) - (eosmiddle->C(6)/eosmiddle->C(0))*(uR_1 - S), eosmiddle->C(4));
						double uRfan = eosmiddle->C(5)*(-eosmiddle->C(0) + eosmiddle->C(7)*uR_1 + S);
						double pRfan = pR_1*pow(eosmiddle->C(5) - (eosmiddle->C(6)/eosmiddle->C(0))*(uR_1 - S), eosmiddle->C(3));
						vector WRfan(dRfan, uRfan, pRfan);
						W.row(i) = WRfan;
					}
				}
			}

		}

	}

	//-------------------------------
	//	output
	//-------------------------------
	std::ofstream outfile;
	outfile.open("dataexact.txt");

	for (int i=0; i<newN; i++){
		double xPos = (static_cast<double>(i) + 0.5)*newdx;
		double e;
		if (Test.tstop >= t_shock){
			double S = (xPos - Test.x1)/t_shock;
			if (S <= ustar){
				e = W(i, 2)/(W(i, 0)*(yM1-1));
			}

			else {
				e = W(i, 2)/(W(i, 0)*(yR-1));
			}
		}
		else{
			double S = (xPos - Test.x0)/Test.tstop;
			if (S <= ustar_1){
				e = W(i, 2)/(W(i, 0)*(yL-1));
			}

			else {
				e = W(i, 2)/(W(i, 0)*(yM1-1));
			}
		}

		outfile << xPos << '\t' << W(i, 0) << '\t' << W(i, 1)
				<< '\t' << W(i, 2) << '\t' << e << std::endl;
	}
	outfile.close();
	std::cout << "done: exact (2 discontinuities)" << std::endl;
}

///FOR STIFFENED GAS////////////////////////////

double GhostFluidMethods::compute_star_pressure_SG(EOS* eosleft, EOS* eosright){
//------------------------------------------
	double TOL = 1e-6;
	auto relative_pressure_change = [](double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
		double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
		return CHA;
	};	
//-------------------------------------------
	//check positivity condition
	try {
		check_pressure_pos_condition(eosleft, eosright);
	} 
	catch(const char* c){
		std::cout << c << std::endl;
		std::cout << "vacuum generated, terminating program" << std::endl;
	}

	StiffenedGas* SGl = dynamic_cast<StiffenedGas*>(eosleft);
	StiffenedGas* SGr = dynamic_cast<StiffenedGas*>(eosright);

	if(SGl == NULL || SGr == NULL) { 
		throw "StiffenedGas exact solver called on non SG EOS";
	}

	double pPV, p0;
	
	//An approximation for p, p0 is required for the initial guess.
	//A poor choice of p0 results in the need for large number of iterations to achieve convergence
	double soundspeed1 = SGl->C(0);
	double soundspeed2 = SGr->C(0);
	double PstiffL = SGl->C(13) + SGl->Pref;
	double PstiffR = SGr->C(13) + SGr->Pref;


	pPV = 0.5*(PstiffL + PstiffR) - (1./8.)*(SGr->C(12) - SGl->C(12))*(SGl->C(11) + SGr->C(11))*(soundspeed1 + soundspeed2);
	pPV = fmax(TOL, pPV);
	double Pmax = fmax(PstiffL, PstiffR);
	double Pmin = fmin(PstiffL, PstiffR);
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
		double y1 = SGl->C(10); double y2 = SGr->C(10);
		if (y1 == y2){
			double pTR = pow((soundspeed1 + soundspeed2 - 0.5*(y1-1)*(SGr->C(12) - SGl->C(12)))/((soundspeed1/pow(PstiffL, SGl->C(1))) + (soundspeed2/pow(PstiffR, SGr->C(1)))), SGl->C(3));
			p0 = pTR;
		}
		else {
			throw "Exact solver attempted on two rarefraction waves with an EOS jump.";
		}
	}
	//If the pressure difference is large or if the guess is larger than Pmax
	else {
		//Select the two shock initial guess using pPV as estimate
		double Q1 = sqrt(SGl->C(8)/(pPV + SGl->C(9)));
		double Q2 = sqrt(SGr->C(8)/(pPV + SGr->C(9)));
		double pTS = (Q1*PstiffL + Q2*PstiffR - (SGr->C(12) - SGr->C(12)))/(Q1 + Q2);
		p0 = pTS;
	}

	/*double pPV, p0;
	double soundspeed1 = eosleft->C(0);
	double soundspeed2 = eosright->C(0);
	pPV = 0.5*(eosleft->C(13) + eosright->C(13)) - (1./8.)*(eosright->C(12) - eosleft->C(12))*(eosleft->C(11) + eosright->C(11))*(soundspeed1 + soundspeed2);
	pPV = fmax(TOL, pPV);
	double Pmax = fmax(eosleft->C(13), eosright->C(13));
	double Pmin = fmin(eosleft->C(13), eosright->C(13));
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
		double y1 = eosleft->C(10); double y2 = eosright->C(10);
		if (y1 == y2){
			double pTR = pow((soundspeed1 + soundspeed2 - 0.5*(y1-1)*(eosright->C(12) - eosleft->C(12)))/((soundspeed1/pow(eosleft->C(13), eosleft->C(1))) + (soundspeed2/pow(eosright->C(13), eosright->C(1)))), eosleft->C(3));
			p0 = pTR;
		}
		else {
			throw "Exact solver attempted on two rarefraction waves with an EOS jump.";
		}
	}
	//If the pressure difference is large or if the guess is larger than Pmax
	else {
		//Select the two shock initial guess using pPV as estimate
		double Q1 = sqrt(eosleft->C(8)/(pPV + eosleft->C(9)));
		double Q2 = sqrt(eosright->C(8)/(pPV + eosright->C(9)));
		double pTS = (Q1*eosleft->C(13) + Q2*eosright->C(13) - (eosright->C(12) - eosright->C(12)))/(Q1 + Q2);
		p0 = pTS;
	}*/

	std::cout << "initial pressure = "<< p0 << std::endl;
	double Pk; Pk = 1e+07;//p0; //First guess
	double Pk_1;
	double CHA = 0;

	int count = 0;
	//std::cout << p0 << std::endl;
	do{
		Pk_1 = newton_raphson(Pk, eosleft, eosright);
		std::cout << "Pk =" << Pk  << '\t' << Pk_1 << std::endl;
		CHA = relative_pressure_change(Pk_1, Pk);
		Pk = Pk_1; //Set the iterate as the new guess
		//std::cout << CHA << '\t' << Pk << std::endl;
		count += 1;

		if (CHA < TOL) break;
		if (count == 20) std::cout << "Warning, maximum iterations reached for Newton's method" << std::endl;
	}while(count < 20);

	std::cout << count <<std::endl;
	return Pk;
}

void GhostFluidMethods::exact_solver_SG(gfmTests Test, EOS* eosleft, EOS* eosright){
	//If there are 2 rarefraction waves, the exact solution no longer holds.

	//setting material EOS parameter
	//eosleft->y = Test.yL;
	//eosright->y = Test.yR;

	StiffenedGas* SGl = dynamic_cast<StiffenedGas*>(eosleft);
	StiffenedGas* SGr = dynamic_cast<StiffenedGas*>(eosright);

	if(SGl == NULL || SGr == NULL) { 
		throw "StiffenedGas exact solver called on non SG EOS";
	}

	vector WL = Test.initialL;	
	vector WR = Test.initialR;

	SGl->y_constants(WL);
	eosright->y_constants(WR);

	double yL = SGl->C(10);
	double yR = SGr->C(10);

	//-----------------------------------------------------------
	//Let the exact solution have a resolution of 1000 grid points
	int newN = 1000;
	double newdx = Test.L/static_cast<double>(newN);
	matrix W(newN, 3);
	//-----------------------------------------------------------

	if (Test.number_of_materials == 2){

		//-----------------Sampling------------------
		//Solution of wavestructure within the domain of interest
		double pstar = compute_star_pressure_SG(eosleft, eosright);
		double ustar = compute_star_velocity(pstar, eosleft, eosright);
		double dshockL = compute_shock_density(pstar, eosleft);
		double dshockR = compute_shock_density(pstar, eosright);
		double drareL = compute_rarefraction_density(pstar, eosleft);
		double drareR = compute_rarefraction_density(pstar, eosright);

		//primitive variables of left and right states
		double dL = WL(0); double dR = WR(0);
		double uL = WL(1); double uR = WR(1);
		double pL = WL(2); double pR = WR(2);

		//Shock Speeds
		double SL = uL - SGl->C(0)*sqrt((SGl->C(2)*(pstar/pL) + SGl->C(1))); 
		double SR = uR + SGr->C(0)*sqrt((SGr->C(2)*(pstar/pR) + SGr->C(1)));

		//Rarefraction wave speeds
		double cstarL = SGl->C(0)*pow(pstar/pL, SGl->C(1));
		double cstarR = SGr->C(0)*pow(pstar/pR, SGr->C(1));

		double SHL = uL - SGl->C(0); double SHR = uR + SGr->C(0); //Head of fan
		double STL = ustar - cstarL; double STR = ustar + cstarR; //tail of fan
		
		//std::cout << drareL << '\t' << drareR << std::endl;
		//std::cout << "Shockspeeds = ";
		//std::cout << SL << '\t' << SR << std::endl;
		//std::cout << SHL << '\t' << SHR << std::endl;
		//std::cout << STL << '\t' << STR << std::endl; 

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
						double dLfan = dL*pow(SGl->C(5) + (SGl->C(6)/SGl->C(0))*(uL - S), SGl->C(4));
						double uLfan = SGl->C(5)*(SGl->C(0) + SGl->C(7)*uL + S);
						double pLfan = pL*pow(SGl->C(5) + (SGl->C(6)/SGl->C(0))*(uL - S), SGl->C(3));
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
						double dRfan = dR*pow(SGr->C(5) - (SGr->C(6)/SGr->C(0))*(uR - S), SGr->C(4));
						double uRfan = SGr->C(5)*(-SGr->C(0) + SGr->C(7)*uR + S);
						double pRfan = pR*pow(SGr->C(5) - (SGr->C(6)/SGr->C(0))*(uR - S), SGr->C(3));
						vector WRfan(dRfan, uRfan, pRfan);
						W.row(i) = WRfan;
					}
				}
			}			
		}

		std::cout << "Pstar = " << pstar << '\t' << "ustar = " << ustar << std::endl;

		//-------------------------------
		//	output
		//-------------------------------
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

	else{
		throw "Exact solver is only avaliable for a single material interface.";
	}

}











