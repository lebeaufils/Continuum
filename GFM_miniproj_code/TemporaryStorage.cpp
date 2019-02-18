double GhostFluidMethods::compute_star_pressure_SG(StiffenedGas SG1, StiffenedGas SG2){
//------------------------------------------
	double TOL = 1e-6;
	auto relative_pressure_change = [](double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
		double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
		return CHA;
	};	
//-------------------------------------------

	EOS* eosleft = &SG1;
	EOS* eosright = &SG2;

	//check positivity condition
	try {
		check_pressure_pos_condition(eosleft, eosright);
	} 
	catch(const char* c){
		std::cout << c << std::endl;
		std::cout << "vacuum generated, terminating program" << std::endl;
	}

	/*double pPV, p0;
	double soundspeed1 = eosleft->C(0);
	double soundspeed2 = eosright->C(0);
	pPV = 0.5*(eosleft->C(13) + eosright->C(13)) - (1./8.)*(eosright->C(12) - eosleft->C(12))*(eosleft->C(11) + eosright->C(11))*(soundspeed1 + soundspeed2);
	pPV = fmax(TOL, pPV);
	double Pmax = fmax(eosleft->C(13), eosright->C(13));
	double Pmin = fmin(eosleft->C(13), eosright->C(13));
	double Quser = 2.0;

	//std::cout << "pPV pressure = "<< pPV << std::endl;

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

		//std::cout << "pTS = " << pTS << std::endl;
	}*/
///////////////////////VERY UNSTABLE NEWTON METHOD, FAILS IF INITIAL PRESSURE GUESS IS > 4E-7////////////////
	//lets cheat by using the geometric mean
	double p0 = sqrt(SG1.C(13)*SG2.C(13));

	//std::cout << "initial pressure = "<< p0 << std::endl;
	double Pk; Pk = p0; //First guess 1e6;//
	double Pk_1;
	double CHA = 0;

	int count = 0;
	//std::cout << p0 << std::endl;
	do{
		Pk_1 = newton_raphson(Pk, eosleft, eosright);
		//std::cout << "Pk =" << Pk  << '\t' << Pk_1 << std::endl;
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

void GhostFluidMethods::exact_solver_SG(gfmTests Test){
	//If there are 2 rarefraction waves, the exact solution no longer holds.
	StiffenedGas SGl;
	StiffenedGas SGr;

	EOS* eosleft = &SGl;
	EOS* eosright = &SGr;

	//setting material EOS parameter
	eosleft->y = Test.yL;
	eosright->y = Test.yR;

	SGl.Pref = Test.Pref1;
	SGr.Pref = Test.Pref2;

	vector WL = Test.initialL;	
	vector WR = Test.initialR;

	eosleft->y_constants(WL);
	eosright->y_constants(WR);

	double yL = eosleft->C(10);
	double yR = eosright->C(10);

	//-----------------------------------------------------------
	//Let the exact solution have a resolution of 1000 grid points
	int newN = 1000;//1000;
	double newdx = Test.L/static_cast<double>(newN);
	matrix W(newN, 3);
	//-----------------------------------------------------------

	if (Test.number_of_materials == 2){

		//-----------------Sampling------------------
		//Solution of wavestructure within the domain of interest
		double pstar = compute_star_pressure_SG(SGl, SGr);
		double ustar = compute_star_velocity(pstar, eosleft, eosright);
		double dshockL = compute_shock_density_SG(pstar, SGl);
		double dshockR = compute_shock_density_SG(pstar, SGr);
		double drareL = compute_rarefraction_density_SG(pstar, SGl);
		double drareR = compute_rarefraction_density_SG(pstar, SGr);

		//primitive variables of left and right states
		double dL = WL(0); double dR = WR(0);
		double uL = WL(1); double uR = WR(1);
		double pL = WL(2); double pR = WR(2);

		//Shock Speeds
		double SL = uL - SGl.C(0)*sqrt((SGl.C(2)*((pstar + SGl.Pref)/(pL + SGl.Pref)) + SGl.C(1))); 
		double SR = uR + SGr.C(0)*sqrt((SGr.C(2)*((pstar + SGr.Pref)/(pR + SGr.Pref)) + SGr.C(1)));

		//Rarefraction wave speeds
		double cstarL = SGl.C(0)*pow((pstar + SGl.Pref)/(pL + SGl.Pref), SGl.C(1));
		double cstarR = SGr.C(0)*pow((pstar + SGr.Pref)/(pR + SGr.Pref), SGr.C(1));

		double SHL = uL - SGl.C(0); double SHR = uR + SGr.C(0); //Head of fan
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
						//std::cout << i << '\t' << "WL" << std::endl;
						W.row(i) = WL;
					}
					//Right of Shock Wave (shocked material, star state)
					else {
						//std::cout << i << '\t' << "WLstar_shock" << std::endl;
						W(i, 0) = dshockL;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
					}
				}
				//Left Rarefraction fan
				else {
					//Left of fastest rarefraction wave (unaffected by rarefraction)
					if (S < SHL){
						//std::cout << i << '\t' << "WL" << std::endl;
						W.row(i) = WL;
					}
					//Out of the rarefraction fan, within the left star state
					else if (S > STL){
						//std::cout << i << '\t' << "WLstar_rare" << std::endl;
						W(i, 0) = drareL;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
					}
					//Within the rarefraction fan
					else {
						//std::cout << i << '\t' << "WLfan" << std::endl;
						double dLfan = dL*pow(SGl.C(5) + (SGl.C(6)/SGl.C(0))*(uL - S), SGl.C(4));
						double uLfan = SGl.C(5)*(SGl.C(0) + SGl.C(7)*uL + S);
						double pLfan = (pL + SGl.Pref)*pow(SGl.C(5) + (SGl.C(6)/SGl.C(0))*(uL - S), SGl.C(3)) - SGl.Pref;
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
						//std::cout << i << '\t' << "WR" << std::endl;
						W.row(i) = WR;
					}
					//Left of Shock Wave (shocked material, star state)
					else {
						//std::cout << i << '\t' << "WRstar_shock" << std::endl;
						W(i, 0) = dshockR;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
					}
				}
				//Right Rarefraction fan
				else {
					//Right of fastest rarefraction wave (unaffected by rarefraction)
					if (S > SHR){
						//std::cout << i << '\t' << "WR" << std::endl;
						W.row(i) = WR;
					}
					//Out of the rarefraction fan, within the right star state
					else if (S < STR){
						//std::cout << i << '\t' << "WRstar_rare" << std::endl;
						W(i, 0) = drareR;
						W(i, 1) = ustar;
						W(i, 2) = pstar;
					}
					//Within the rarefraction fan
					else {
						//std::cout << i << '\t' << "WRfan" << std::endl;
						double dRfan = dR*pow(SGr.C(5) - (SGr.C(6)/SGr.C(0))*(uR - S), SGr.C(4));
						double uRfan = SGr.C(5)*(-SGr.C(0) + SGr.C(7)*uR + S);
						double pRfan = pR*pow(SGr.C(5) - (SGr.C(6)/SGr.C(0))*(uR - S), SGr.C(3));
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
				e = (W(i, 2) + SGl.y*SGl.Pref)/(W(i, 0)*(yL-1));
			}

			else {
				e = (W(i, 2) + SGr.y*SGr.Pref)/(W(i, 0)*(yR-1));
			}

			outfile << xPos << '\t' << W(i, 0) << '\t' << W(i, 1)
					<< '\t' << W(i, 2) << '\t' << e << std::endl;
		}
		outfile.close();
		std::cout << "done: exact" << std::endl;
	}

	else{
		throw "Exact solver for StiffenedGas is only avaliable for a single material interface.";
	}

}

//testing
double GhostFluidMethods::compute_shock_density_SG(double pstar, StiffenedGas eos){
	//turn this into a stored variable so the iteration does not have to be called multiple times
	//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
	double Pratio = (pstar + eos.Pref)/(eos.C(13) + eos.Pref);
	double dshock = eos.C(11)*((eos.C(6) + Pratio)/(eos.C(6)*Pratio + 1));
	return dshock;
}

double GhostFluidMethods::compute_rarefraction_density_SG(double pstar, StiffenedGas eos){
	double drare = eos.C(11)*pow((pstar + eos.Pref)/(eos.C(13) + eos.Pref), 1./eos.C(10));
	return drare;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

void GhostFluidMethods::ghost_boundary_RP_SG(MUSCL* Uleft, StiffenedGas SGl, MUSCL* Uright, StiffenedGas SGr, int i){
//Exact riemannsolver for idealgas and stiffened gas. Ghost point values are the RP star states.
//Consider an interface between i and i+1
	//converting to primitive variables

	double dl = Uleft->U(i-1, 0);
	double dr = Uright->U(i+2, 0);

	double ul = Uleft->U(i-1, 1)/Uleft->U(i-1, 0);
	double ur = Uright->U(i+2, 1)/Uright->U(i+2, 0);

	double Pl = SGl.Pressure(Uleft->U, i-1);
	double Pr = SGr.Pressure(Uright->U, i+2);

	vector WL(dl, ul, Pl);
	vector WR(dr, ur, Pr);

	SGl.y_constants(WL);
	SGr.y_constants(WR);

	EOS* eosleft = &SGl;
	EOS* eosright = &SGr;

	//Solution of wavestructure within the domain of interest
	double pstar = compute_star_pressure_SG(SGl, SGr);
	//std::cout << "pstar = " << pstar << std::endl;
	double ustar = compute_star_velocity(pstar, eosleft, eosright);
	double dshockL = compute_shock_density_SG(pstar, SGl);
	double dshockR = compute_shock_density_SG(pstar, SGr);
	double drareL = compute_rarefraction_density_SG(pstar, SGl);
	double drareR = compute_rarefraction_density_SG(pstar, SGr);

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
	Uleft->U.row(i) = SGl.conservedVar(WLghost);
	Uright->U.row(i+1) = SGr.conservedVar(WRghost);
	for (int j=0; j<4; j++){
		//Defining the states at ghost fluid cells
		Uleft->U.row(i+j+1) = SGl.conservedVar(WLghost);
		Uright->U.row(i-j) = SGr.conservedVar(WRghost);
	}
}


void GhostFluidMethods::initial_conditions_RP_SG(StiffenedGas SGl, StiffenedGas SGr, gfmTests Test){

	vector euler1;
	vector euler2;

	EOS* eos1 = &SGl; EOS* eos2 = &SGr;

	//setting material EOS parameter
	eos1->y = Test.yL; SGl.Pref = Test.Pref1;
	eos2->y = Test.yR; SGr.Pref = Test.Pref2;

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
			ghost_boundary_RP_SG(var1, SGl, var2, SGr, i+1); //matrix 2 contains ghost points
		}
	}

	for (int i=0; i<N; i++){
		//std::cout << i << '\t' << phi(i+1) << '\t' << get_sgn(phi(i+1)) << var1->U.row(i+2) << '\t' << '\t' << var2->U.row(i+2) << std::endl;
	}

	var1->boundary_conditions();
	var2->boundary_conditions();
	boundary_conditions(); //levelsetfunction
}

void GhostFluidMethods::solver_RP_SG(gfmTests Test){
	slopeLimiter a;
	a = var1->getLimiter();

	StiffenedGas SGl;
	StiffenedGas SGr;

	EOS* eos1 = &SGl;
	EOS* eos2 = &SGr;

	initial_conditions_RP_SG(SGl, SGr, Test);

	double t = 0.0;

	do{
		//compute ghost fluid boundaries
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){
				//std::cout << count << " boundary is at i = " << i << std::endl;
				ghost_boundary_RP_SG(var1, SGl, var2, SGr, i+1); //matrix 2 contains ghost points
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
			if (phi(i-1) > 0){
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
		//if (count == 50) t = Test.tstop;

	}while(t<Test.tstop);
	std::cout << count << std::endl;
	output(eos1, eos2);

}
///////////////////////////////////////////////////////////////////






////////////////////////////////////////
	StiffenedGas* SGl = dynamic_cast<StiffenedGas*>(eos1);
	StiffenedGas* SGr = dynamic_cast<StiffenedGas*>(eos2);

	if (SGl == NULL && SGr == NULL){

		//After setting the real fluids, populate the empty cells with ghost values
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){
				//std::cout << "boundary is at i = " << i << std::endl;
				ghost_boundary_RP(var1, eos1, var2, eos2, i+1); //matrix 2 contains ghost points
			}
		}

	}

	else if (SGl == NULL || SGr == NULL){
		throw "Solver attempted on two different equations of states";
	}

	else {
		for (int i=1; i<N+1; i++){
			int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
			//std::cout << phi(i) << '\t' << testsgn << std::endl;
			if (testsgn > -2 && testsgn < 1 && i!=N){
				//std::cout << "boundary is at i = " << i << std::endl;
				ghost_boundary_RP_SG(var1, SGl, var2, SGr, i+1); //matrix 2 contains ghost points
			}
		}
	}



	StiffenedGas* SGl = dynamic_cast<StiffenedGas*>(eos1);
	StiffenedGas* SGr = dynamic_cast<StiffenedGas*>(eos2);

	if (SGl == NULL && SGr == NULL){	
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
//If it is stiffened gas
	else {
		do{
			//compute ghost fluid boundaries
			for (int i=1; i<N+1; i++){
				int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
				//std::cout << phi(i) << '\t' << testsgn << std::endl;
				if (testsgn > -2 && testsgn < 1 && i!=N){
					//std::cout << count << " boundary is at i = " << i << std::endl;
					ghost_boundary_RP_SG(var1, SGl, var2, SGr, i+1); //matrix 2 contains ghost points
				}

				if (count==0) std::cout << var1->U.row(i) << std::endl; 
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
//-----------------------------------------------------------------
//	EXACT
//-----------------------------------------------------------------
void GhostFluidMethods::y_constants(gfmTests Test){
	for (int i=0; i<Test.number_of_materials; i++){
		// The order of materials stored is
		// Left, Right, Middle (1), Middle (2)
		double y;
		vector W;

		if (i == 0) {
			y = Test.yL;
			W = Test.initialL;
		}
		else if (i == 1) {
			y = Test.yR;
			W = Test.initialR;
		}
		else if (i == 2) {
			y = Test.yM1;
			W = Test.initialM1;
		}
		else if (i == 3) {
			y = Test.yM2;
			W = Test.initialM2;
		}
		else {
			throw "Exceeded maximum allowable materials";
		} 

		C(i, 0) = sqrt(y*W(2)/W(0));
		C(i, 1) = (y-1)/(2*y); // y-1 / 2y
		C(i, 2) = (y+1)/(2*y); // y+1 / 2y
		C(i, 3) = (2*y)/(y-1);
		C(i, 4) = 2./(y-1);
		C(i, 5) = 2./(y+1);
		C(i, 6) = (y-1)/(y+1);
		C(i, 7) = (y-1)/2.;
		C(i, 8) = (2./(y+1))/W(0);//Ak
		C(i, 9) = W(2)*((y-1)/(y+1));  //Bk
		C(i, 10) = y;
	}
}

//note: the left material is always listed first
double GhostFluidMethods::f(double P, EOS* eosleft, EOS* eosright){
	double du = (eosright->C(12) - eos1->C(12)); //velocity difference
	return eos1->fk(P) + eos2->fk(P) + du;   
}

double GhostFluidMethods::fprime(double P, vector W1, int mat1, vector W2, int mat2){
	return fprimek(P, W1, mat1) + fprimek(P, W2, mat2);  
}

void GhostFluidMethods::check_pressure_pos_condition(vector W1, int mat1, vector W2, int mat2) { 
	//(Δu)crit ≡ 2aL/γ−1 + 2aR/γ−1 ≤ uR −uL , (4.82) Toro pg 127
	double du_crit = C(mat1, 4)*C(mat1, 0) + C(mat2, 4)*C(mat2, 0);
	double du = W2(1) - W1(1);

	//If the pressure positivity condition is not satisfied, vacuum is generated
	// and the solver fails. Exit the program if this condition is not satisfied.
	if (du_crit <= du){ //This conditions ensures S*L <= S*R
		throw "Pressure positivity condition violated";
	}
}

double GhostFluidMethods::newton_raphson(double Pk, vector W1, int mat1, vector W2, int mat2){
	double Pk_1 = Pk - f(Pk, W1, mat1, W2, mat2)/fprime(Pk, W1, mat1, W2, mat2);
	return Pk_1;
}

double GhostFluidMethods::compute_star_pressure(vector W1, int mat1, vector W2, int mat2){
//------------------------------------------
	double TOL = 1e-6;
	auto relative_pressure_change = [](double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
		double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
		return CHA;
	};	
//-------------------------------------------
	//check positivity condition
	try {
		check_pressure_pos_condition(W1, mat1, W2, mat2);
	} 
	catch(const char* c){
		std::cout << c << std::endl;
		std::cout << "vacuum generated, terminating program" << std::endl;
	}
	double pPV, p0;
	
	//An approximation for p, p0 is required for the initial guess.
	//A poor choice of p0 results in the need for large number of iterations to achieve convergence
	double soundspeed1 = sqrt(C(mat1, 10)*W1(2)/W1(0));
	double soundspeed2 = sqrt(C(mat1, 10)*W2(2)/W2(0));
	pPV = 0.5*(W1(2) + W2(2)) - (1./8.)*(W2(1) - W1(1))*(W1(0) + W2(0))*(soundspeed1 + soundspeed2);
	pPV = fmax(TOL, pPV);
	double Pmax = fmax(W1(2), W2(2));
	double Pmin = fmin(W1(2), W2(2));
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
		double y1 = C(mat1, 10); double y2 = C(mat2, 10);
		if (y1 == y2){
			double pTR = pow((soundspeed1 + soundspeed2 - 0.5*(y1-1)*(W2(1) - W1(1)))/((soundspeed1/pow(W1(2), C(mat1, 1))) + (C(1, 0)/pow(W2(2), C(0, 1)))), C(0, 3));
			p0 = pTR;
		}
		else {
			throw "Exact solver attempted on two rarefraction waves with an EOS jump.";
		}
	}
	//If the pressure difference is large or if the guess is larger than Pmax
	else {
		//Select the two shock initial guess using pPV as estimate
		double Q1 = sqrt(C(mat1, 8)/(pPV + C(mat1, 9)));
		double Q2 = sqrt(C(mat2, 8)/(pPV + C(mat2, 9)));
		double pTS = (Q1*W1(2) + Q2*W2(2) - (W2(1) - W1(1)))/(Q1 + Q2);
		p0 = pTS;
	}

	double Pk; Pk = p0; //First guess
	double Pk_1;
	double CHA = 0;
	double mixTOL;

	int count = 0;
	//std::cout << p0 << std::endl;
	do{
		Pk_1 = newton_raphson(Pk, W1, mat1, W2, mat2);
		CHA = relative_pressure_change(Pk_1, Pk);
		Pk = Pk_1; //Set the iterate as the new guess
		//std::cout << CHA << '\t' << Pk << std::endl;
		count += 1;
		mixTOL = TOL*(1 + Pk);

		if (CHA < mixTOL) break;
		if (count == 20) std::cout << "Warning, maximum iterations reached for Newton's method" << std::endl;
	}while(count < 20);

	//std::cout << count <<std::endl;
	return Pk;
}

double GhostFluidMethods::compute_star_velocity(double pstar, vector W1, int mat1, vector W2, int mat2){
	double ustar = 0.5*(W1(1) + W2(1)) + 0.5*(fk(pstar, W2, mat2) - fk(pstar, W1, mat1));
	return ustar;
}

double GhostFluidMethods::compute_shock_density(double pstar, vector Wk, int mat){
	//turn this into a stored variable so the iteration does not have to be called multiple times
	//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
	double Pratio = pstar/Wk(2);
	double dshock = Wk(0)*((C(mat, 6) + Pratio)/(C(mat, 6)*Pratio + 1));
	return dshock;
}

double GhostFluidMethods::compute_rarefraction_density(double pstar, vector Wk, int mat){
	double drare = Wk(0)*pow(pstar/Wk(2), 1./C(mat, 10));
	return drare;
}


void GhostFluidMethods::exact_solver(gfmTests Test){
	y_constants(Test);
	//If there are 2 rarefraction waves, the exact solution no longer holds.
	double yL = C(0, 10);
	double yR = C(1, 10);

	vector WL = Test.initialL; //[0]
	vector WR = Test.initialR; //[1]

	//-----------------------------------------------------------
	//Let the exact solution have a resolution of 1000 grid points
	int newN = 1000;
	double newdx = Test.L/static_cast<double>(newN);
	matrix W(newN, 3);
	//-----------------------------------------------------------

	if (Test.number_of_materials == 2){

		//-----------------Sampling------------------
		//Solution of wavestructure within the domain of interest
		double pstar = compute_star_pressure(WL, 0, WR, 1);
		double ustar = compute_star_velocity(pstar, WL, 0, WR, 1);
		double dshockL = compute_shock_density(pstar, WL, 0);
		double dshockR = compute_shock_density(pstar, WR, 1);
		double drareL = compute_rarefraction_density(pstar, WL, 0);
		double drareR = compute_rarefraction_density(pstar, WR, 1);

		//primitive variables of left and right states
		double dL = WL(0); double dR = WR(0);
		double uL = WL(1); double uR = WR(1);
		double pL = WL(2); double pR = WR(2);

		//Shock Speeds
		double SL = uL - C(0, 0)*sqrt((C(0, 2)*(pstar/pL) + C(0, 1))); 
		double SR = uR + C(1, 0)*sqrt((C(1, 2)*(pstar/pR) + C(1, 1)));

		//Rarefraction wave speeds
		double cstarL = C(0, 0)*pow(pstar/pL, C(0, 1));
		double cstarR = C(1, 0)*pow(pstar/pR, C(1, 1));

		double SHL = uL - C(0, 0); double SHR = uR + C(1, 0); //Head of fan
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
						double dLfan = dL*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL - S), C(0, 4));
						double uLfan = C(0, 5)*(C(0, 0) + C(0, 7)*uL + S);
						double pLfan = pL*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL - S), C(0, 3));
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
						double dRfan = dR*pow(C(1, 5) - (C(1, 6)/C(1, 0))*(uR - S), C(1, 4));
						double uRfan = C(1, 5)*(-C(1, 0) + C(1, 7)*uR + S);
						double pRfan = pR*pow(C(1, 5) - (C(1, 6)/C(1, 0))*(uR - S), C(1, 3));
						vector WRfan(dRfan, uRfan, pRfan);
						W.row(i) = WRfan;
					}
				}
			}			
		}

		std::cout << "Pstar = " << pstar << '\t' << "ustar = " << ustar << std::endl;

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

	//NOTE: only works on Fedkiw (2002), Test B and C
	else if (Test.number_of_materials == 3){
		//For 3 masterisls, Consider 2 RP between Left and Middle states, as well as Right and Middle states.
		//Assuming that the stop time is before interactions between the generated waves

		double yM1 = C(2, 10);
		vector WM1 = Test.initialM1; //[2]

		//----------------Rename variables--------------------------
		//variablename_1 represents solution to the leftmost RP
		//----------------------------------------------------------
		//Solution of wavestructure within the domain of interest between WL and WM1
		//This is shared with the solver for 4 initial conditions
		double pstar_1 = compute_star_pressure(WL, 0, WM1, 2);
		double ustar_1 = compute_star_velocity(pstar_1, WL, 0, WM1, 2);
		double dshockL_1 = compute_shock_density(pstar_1, WL, 0);
		double dshockR_1 = compute_shock_density(pstar_1, WM1, 2);
		double drareL_1 = compute_rarefraction_density(pstar_1, WL, 0);
		double drareR_1 = compute_rarefraction_density(pstar_1, WM1, 2);

		//primitive variables of left and right states
		double dL_1 = WL(0); double dR_1 = WM1(0);
		double uL_1 = WL(1); double uR_1 = WM1(1);
		double pL_1 = WL(2); double pR_1 = WM1(2);

		//Shock Speeds
		double SL_1 = uL_1 - C(0, 0)*sqrt((C(0, 2)*(pstar_1/pL_1) + C(0, 1))); 
		double SR_1 = uR_1 + C(2, 0)*sqrt((C(2, 2)*(pstar_1/pR_1) + C(2, 1)));

		//Rarefraction wave speeds
		double cstarL_1 = C(0, 0)*pow(pstar_1/pL_1, C(0, 1));
		double cstarR_1 = C(2, 0)*pow(pstar_1/pR_1, C(2, 1));

		double SHL_1 = uL_1 - C(0, 0); double SHR_1 = uR_1 + C(2, 0); //Head of fan
		double STL_1 = ustar_1 - cstarL_1; double STR_1 = ustar_1 + cstarR_1; //tail of fan

	//Tracking the right moving shock wave. For this problem, we know that the right moving wave of the
		//first riemann problem is a shock wave.
		//time when shock wave reaches the discontinuity
		double t_shock = Test.tstop - (Test.x1 - Test.x0)/SR_1;
		vector WstarL(dshockR_1, ustar_1, pstar_1);	

		//Solution of wavestructure within the domain of interest
		double pstar = compute_star_pressure(WstarL, 2, WR, 1); //WstarL using y of the middle state
		double ustar = compute_star_velocity(pstar, WstarL, 2, WR, 1);
		double dshockL = compute_shock_density(pstar, WstarL, 2);
		double dshockR = compute_shock_density(pstar, WR, 1);
		double drareL = compute_rarefraction_density(pstar, WstarL, 2);
		double drareR = compute_rarefraction_density(pstar, WR, 1);

		//primitive variables of left and right states
		double dL = WstarL(0); double dR = WR(0);
		double uL = WstarL(1); double uR = WR(1);
		double pL = WstarL(2); double pR = WR(2);

		double soundspeedL = sqrt(C(2, 10)*WstarL(2)/WstarL(0));
		double soundspeedR = C(1, 0);

		//Shock Speeds
		double SL = uL - soundspeedL*sqrt((C(2, 2)*(pstar/pL) + C(2, 1))); 
		double SR = uR + soundspeedR*sqrt((C(1, 2)*(pstar/pR) + C(1, 1)));

		//Rarefraction wave speeds
		double cstarL = soundspeedL*pow(pstar/pL, C(2, 1));
		double cstarR = soundspeedR*pow(pstar/pR, C(1, 1));

		double SHL = uL - soundspeedL; double SHR = uR + soundspeedR; //Head of fan
		double STL = ustar - cstarL; double STR = ustar + cstarR; //tail of fan
	

		//-----------------------------------------------------------
		//	Samp[ling]
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
							double dLfan = dL*pow(C(2, 5) + (C(2, 6)/soundspeedL)*(uL - S), C(2, 4));
							double uLfan = C(2, 5)*(soundspeedL + C(2, 7)*uL + S);
							double pLfan = pL*pow(C(2, 5) + (C(2, 6)/soundspeedL)*(uL - S), C(2, 3));
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
							double dRfan = dR*pow(C(2, 5) - (C(2, 6)/soundspeedR)*(uR - S), C(2, 4));
							double uRfan = C(2, 5)*(-soundspeedR + C(2, 7)*uR + S);
							double pRfan = pR*pow(C(2, 5) - (C(2, 6)/soundspeedR)*(uR - S), C(2, 3));
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
							double dLfan = dL_1*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL_1 - S), C(0, 4));
							double uLfan = C(0, 5)*(C(0, 0) + C(0, 7)*uL_1 + S);
							double pLfan = pL_1*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL_1 - S), C(0, 3));
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
							double dRfan = dR_1*pow(C(2, 5) - (C(2, 6)/C(2, 0))*(uR_1 - S), C(2, 4));
							double uRfan = C(2, 5)*(-C(2, 0) + C(2, 7)*uR_1 + S);
							double pRfan = pR_1*pow(C(2, 5) - (C(2, 6)/C(2, 0))*(uR_1 - S), C(2, 3));
							vector WRfan(dRfan, uRfan, pRfan);
							W.row(i) = WRfan;
						}
					}
				}

			}

		}

		/*-------------------------------
			output
		-------------------------------*/
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
	
}