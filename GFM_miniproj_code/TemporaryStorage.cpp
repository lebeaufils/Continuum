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