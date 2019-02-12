/*void GhostFluidMethods::exact_solver_old(gfmTests Test){
	y_constants(Test);
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
	//double C(0, 4) = 2./(yL-1);
	double CONST5L = 2./(yL+1);
	double CONST6L = (yL-1)/(yL+1);
	double CONST7L = (yL-1)/2.;
	double CONST8L = yL-1;

	double CONST1R = (yR-1)/(2*yR); // y-1 / 2y
	double CONST2R = (yR+1)/(2*yR); // y+1 / 2y
	double CONST3R = (2*yR)/(yR-1);
	double CONST4R = 2./(yR-1);
	//double C(1, 5) = 2./(yR+1);
	//double C(1, 6) = (yR-1)/(yR+1);
	double CONST7R = (yR-1)/2.;

	double TOL = 1e-6;

	//data-dependent constants
	double AL = C(0, 5)/WL(0); 	double AR = C(1, 5)/WR(0);
	double BL = WL(2)*C(0, 6); 	double BR = WR(2)*C(1, 6);
	//-------------------------------------------

	//generic functions for 2 or more materials
	//function of P, f(P, WL, WR)
	auto fL = [cL, AL, BL, WL, this, CONST1L](double P){
		double QL = sqrt(AL/(P + BL));
		double FL;

		if (P > WL(2)){ //Shock
			FL = (P - WL(2))*QL;
		}

		else { //Rarefraction
			FL = C(0, 4)*cL*(pow(P/WL(2), CONST1L) - 1);
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

	auto relative_pressure_change = [](double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
		double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
		return CHA;
	};	

	//-------------------------------------------
	if (Test.number_of_materials == 2){
		auto check_pressure_pos_condition = [cL, cR, WL, WR, this, CONST4R]() { 
			//(Δu)crit ≡ 2aL/γ−1 + 2aR/γ−1 ≤ uR −uL , (4.82) Toro pg 127
			double du_crit = C(0, 4)*cL + CONST4R*cR;
			double du = WR(1) - WL(1);

			//If the pressure positivity condition is not satisfied, vacuum is generated
			// and the solver fails. Exit the program if this condition is not satisfied.
			if (du_crit <= du){ //This conditions ensures S*L <= S*R
				throw "Pressure positivity condition violated";
			}
		};

		auto f = [WL, WR, fR, fL](double P){
			double du = (WR(1) - WL(1)); //velocity difference
			return fL(P) + fR(P) + du;   
		};


		auto fprime = [fprimeL, fprimeR](double P){
			return fprimeL(P) + fprimeR(P);  
		};

		auto newton_raphson = [f, fprime](double Pk){
			double Pk_1 = Pk - f(Pk)/fprime(Pk);
			return Pk_1;
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
				std::cout << CHA << '\t' << Pk << std::endl;
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

		auto compute_shock_density_R = [WR, this](double pstar){
			//turn this into a stored variable so the iteration does not have to be called multiple times
			//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
			double Pratio = pstar/WR(2);
			double dshock_R = WR(0)*((C(1, 6) + Pratio)/(C(1, 6)*Pratio + 1));
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

	//-----------------------------------------------------------
	//	Sampling
	//-----------------------------------------------------------
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
		
		//std::cout << drareL << '\t' << drareR << std::endl;
		//std::cout << SHL << '\t' << SHR << std::endl;
		//std::cout << STL << '\t' << STR<< std::endl; 

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
						double dLfan = dL*pow(CONST5L + (CONST6L/cL)*(uL - S), C(0, 4));
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
						double dRfan = dR*pow(C(1, 5) - (C(1, 6)/cR)*(uR - S), CONST4R);
						double uRfan = C(1, 5)*(-cR + CONST7R*uR + S);
						double pRfan = pR*pow(C(1, 5) - (C(1, 6)/cR)*(uR - S), CONST3R);
						vector WRfan(dRfan, uRfan, pRfan);
						W.row(i) = WRfan;
					}
				}
			}			
		}

		std::cout << compute_star_pressure() << std::endl;

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
	/*
	else if (Test.number_of_materials >= 3){
		//For 3 masterisls, Consider 2 RP between Left and Middle states, as well as Right and Middle states.
		//Assuming that the stop time is before interactions between the generated waves
		double yM = Test.y3;

		vector WM = Test.initialM1;

		//-------------------------------------------
		//constants
		double cM = sqrt(yL*WM(2)/WM(0));

		double CONST1M = (yM-1)/(2*yM); // y-1 / 2y
		double CONST2M = (yM+1)/(2*yM); // y+1 / 2y
		double CONST3M = (2*yM)/(yM-1);
		double CONST4M = 2./(yM-1);
		double CONST5M = 2./(yM+1);
		double CONST6M = (yM-1)/(yM+1);
		double CONST7M = (yM-1)/2.;
		double CONST8M = yM-1;

		double TOL = 1e-6;

		//data-dependent constants
		double AM = CONST5M/WM(0);
		double BM = WM(2)*CONST6M; 
		//-------------------------------------------
		//function of P, f(P, WL, WR)
		auto fM = [cM, AM, BM, WM, CONST4M, CONST1M](double P){
			double QM = sqrt(AM/(P + BM));
			double FM;

			if (P > WM(2)){ //Shock
				FM = (P - WM(2))*QM;
			}

			else { //Rarefraction
				FM = CONST4M*cM*(pow(P/WM(2), CONST1M) - 1);
			}

			return FM;
		};

		auto fprimeM = [cM, AM, BM, WM, CONST2M](double P){
			double QM = sqrt(AM/(P + BM)); 
			double FMprime;
			if (P > WM(2)){ //Shock
				FMprime = QM*(1 - (P - WM(2))/(2.*(P + BM))); 
			}

			else { //Rarefraction
				FMprime = (1./(WM(0)*cM))*pow(P/WM(2), -CONST2M);
			}
			return FMprime;
		};

		//Riemann Problem Between Left and Middle states
		auto f_LM = [WL, WM, fL, fM](double P){
			double du = (WM(1) - WL(1)); //velocity difference
			return fL(P) + fM(P) + du;   
		};


		auto fprime_LM = [fprimeL, fprimeM](double P){
			return fprimeL(P) + fprimeM(P);  
		};

		auto newton_raphson_LM = [f_LM, fprime_LM](double Pk){
			double Pk_1 = Pk - f_LM(Pk)/fprime_LM(Pk);
			return Pk_1;
		};

		//Riemann Problem Between Middle and Right states
		auto f_MR = [WM, WR, fM, fR](double P){
			double du = (WR(1) - WM(1)); //velocity difference
			return fM(P) + fR(P) + du;   
		};


		auto fprime_MR = [fprimeM, fprimeR](double P){
			return fprimeM(P) + fprimeR(P);  
		};

		auto newton_raphson_MR = [f_MR, fprime_MR](double Pk){
			double Pk_1 = Pk - f_MR(Pk)/fprime_MR(Pk);
			return Pk_1;
		};

		if (Test.number_of_materials == 3){

		}

		//else if (Test.number_of_materials == 4){}

	}
}*/

				//std::cout << "distance from x0 = "<< (xPos - newx0) << '\t' << "distance from x1 = " << (xPos - newx1)<< std::endl;
				//std::cout << i << '\t' << xPos << std::endl;
				if (xPos <= Pos_R_acoustic){
					double S = (xPos - Test.x0)/Test.tstop;
					//Sampled point lies in the solution of the Left RP
					//std::cout << i << '\t'<< S << std::endl;

					//Left side of Contact wave
					if (S <= ustar_1){
						//Left Shock wave
						if (pstar_1 > pL_1){
							//Left of Shock Wave (unaffected by shock)
							if (S < SL_1){
								W.row(i) = WL;
								std::cout << "Shock WL" << '\t' << i << std::endl;
							}
							//Right of Shock Wave (shocked material, star state)
							else {
								W(i, 0) = dshockL_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;
								std::cout << "Shock WL star" << '\t' << i << std::endl;
							}
						}
						//Left Rarefraction fan
						else {
							//Left of fastest rarefraction wave (unaffected by rarefraction)
							if (S < SHL_1){
								W.row(i) = WL;
								std::cout << "Rarefraction WL" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the left star state
							else if (S > STL_1){
								W(i, 0) = drareL_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;
								std::cout << "Rarefraction WL star" << '\t' << i << std::endl;
							}
							//Within the rarefraction fan
							else {
								double dLfan = dL_1*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL_1 - S), C(0, 4));
								double uLfan = C(0, 5)*(C(0, 0) + C(0, 7)*uL_1 + S);
								double pLfan = pL_1*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL_1 - S), C(0, 3));
								vector WLfan(dLfan, uLfan, pLfan);
								W.row(i) = WLfan;
								std::cout << "Rarefraction WL fan" << '\t' << i << std::endl;
							}
						}
					}
					//right side of contact wave
					else {
						//Right Shock wave
						if (pstar_1 > pR_1){
							//Right of Shock Wave (unaffected by shock)
							if (S > SR_1){
								std::cout << "Shock WM" << '\t' << i << std::endl;
								W.row(i) = WM1;

							}
							//Left of Shock Wave (shocked material, star state)
							else {
								std::cout << "Shock WM Star" << '\t' << i << std::endl;
								W(i, 0) = dshockR_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;

							}
						}
						//Right Rarefraction fan
						else {
							//Right of fastest rarefraction wave (unaffected by rarefraction)
							if (S > SHR_1){
								std::cout << "rarefraction WM" << '\t' << i << std::endl;
								W.row(i) = WM1;
							}
							//Out of the rarefraction fan, within the right star state
							else if (S < STR_1){
								std::cout << "rarefraction WM Star" << '\t' << i << std::endl;
								W(i, 0) = drareR_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								std::cout << "rarefraction WM fan" << '\t' << i << std::endl;
								double dRfan = dR_1*pow(C(2, 5) - (C(2, 6)/C(2, 0))*(uR_1 - S), C(2, 4));
								double uRfan = C(2, 5)*(-C(2, 0) + C(2, 7)*uR_1 + S);
								double pRfan = pR_1*pow(C(2, 5) - (C(2, 6)/C(2, 0))*(uR_1 - S), C(2, 3));
								vector WRfan(dRfan, uRfan, pRfan);
								W.row(i) = WRfan;
							}
						}
					}			
				}

				else if (xPos >= Pos_L_acoustic){ 
					//------------------------------------------------
					// sampled point lies in solution to Right RP
					//------------------------------------------------
					double S = (xPos - Test.x1)/Test.tstop;
					//Left side of Contact wave
					if (S <= ustar_2){
						//Left Shock wave
						if (pstar_2 > pL_2){
							//Left of Shock Wave (unaffected by shock)
							if (S < SL_2){
								W.row(i) = WM1;
								std::cout << "right Shock WM" << '\t' << i << std::endl;
							}
							//Right of Shock Wave (shocked material, star state)
							else {
								W(i, 0) = dshockL_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "right Shock WM Star" << '\t' << i << std::endl;
							}
						}
						//Left Rarefraction fan
						else {
							//Left of fastest rarefraction wave (unaffected by rarefraction)
							if (S < SHL_2){
								W.row(i) = WM1;
								std::cout << "right rarefraction WM" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the left star state
							else if (S > STL_2){
								W(i, 0) = drareL_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "right rarefraction WM Star" << '\t' << i << std::endl;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								double dLfan = dL_2*pow(C(2, 5) + (C(2, 6)/C(2, 0))*(uL_2 - S), C(2, 4));
								double uLfan = C(2, 5)*(C(2, 0) + C(2, 7)*uL_2 + S);
								double pLfan = pL_2*pow(C(2, 5) + (C(2, 6)/C(2, 0))*(uL_2 - S), C(2, 3));
								vector WLfan(dLfan, uLfan, pLfan);
								W.row(i) = WLfan;
								std::cout << "right rarefraction WM fan" << '\t' << i << std::endl;
							}
						}
					}
					//right side of contact wave
					else {
						//Right Shock wave
						if (pstar_2 > pR_2){
							//Right of Shock Wave (unaffected by shock)
							if (S > SR_2){
								//std::cout << "Shock WR" << '\t' << i << std::endl;
								W.row(i) = WR;
								std::cout << "Shock WR" << '\t' << i << std::endl;
							}
							//Left of Shock Wave (shocked material, star state)
							else {
								//std::cout << "Shock Star" << '\t' << i << std::endl;
								W(i, 0) = dshockR_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "Shock WR Star" << '\t' << i << std::endl;
							}
						}
						//Right Rarefraction fan
						else {
							//Right of fastest rarefraction wave (unaffected by rarefraction)
							if (S > SHR_2){
								//std::cout << "rarefraction WR" << '\t' << i << std::endl;
								W.row(i) = WR;
								std::cout << "rarefraction WR" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the right star state
							else if (S < STR_2){
								//std::cout << "rarefraction Star" << '\t' << i << std::endl;
								W(i, 0) = drareR_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "rarefraction WR Star" << '\t' << i << std::endl;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								//std::cout << "rarefraction fan" << '\t' << i << std::endl;
								double dRfan = dR_2*pow(C(1, 5) - (C(1, 6)/C(1, 0))*(uR_2 - S), C(1, 4));
								double uRfan = C(1, 5)*(-C(1, 0) + C(1, 7)*uR_2 + S);
								double pRfan = pR_2*pow(C(1, 5) - (C(1, 6)/C(1, 0))*(uR_2 - S), C(1, 3));
								vector WRfan(dRfan, uRfan, pRfan);
								W.row(i) = WRfan;
								std::cout << "rarefraction WR fan" << '\t' << i << std::endl;
							}
						}
					}					
				}

				else{ //in between the 2 starstates
	
					double S = (xPos - newx1)/Test.tstop;
					//Left side of Contact wave
					if (S <= ustar){
						//Left Shock wave
						if (pstar > pL){
							//Left of Shock Wave (unaffected by shock)
							if (S < SL){
								W.row(i) = WstarL;
								std::cout << "Shock WstarL" << '\t' << i << std::endl;
							}
							//Right of Shock Wave (shocked material, star state)
							else {
								W(i, 0) = dshockL;
								W(i, 1) = ustar;
								W(i, 2) = pstar;
								std::cout << "Shock new WstarL" << '\t' << i << std::endl;
							}
						}
						//Left Rarefraction fan
						else {
							//Left of fastest rarefraction wave (unaffected by rarefraction)
							if (S < SHL){
								W.row(i) = WstarL;
								std::cout << "Rarefraction WstarL" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the left star state
							else if (S > STL){
								W(i, 0) = drareL;
								W(i, 1) = ustar;
								W(i, 2) = pstar;
								std::cout << "Rarefraction new WstarL" << '\t' << i << std::endl;
							}
							//Within the rarefraction fan
							else {
								double dLfan = dL*pow(C(2, 5) + (C(2, 6)/soundspeedL)*(uL - S), C(2, 4));
								double uLfan = C(2, 5)*(soundspeedL + C(2, 7)*uL + S);
								double pLfan = pL*pow(C(2, 5) + (C(2, 6)/soundspeedL)*(uL - S), C(2, 3));
								vector WLfan(dLfan, uLfan, pLfan);
								W.row(i) = WLfan;
								std::cout << "Rarefraction new WfanL" << '\t' << i << std::endl;
							}
						}
					}
					//right side of contact wave
					else {
						//Right Shock wave
						if (pstar > pR){
							//Right of Shock Wave (unaffected by shock)
							if (S > SR){
								std::cout << "Shock WstarR" << '\t' << i << std::endl;
								W.row(i) = WstarR;

							}
							//Left of Shock Wave (shocked material, star state)
							else {
								std::cout << "Shock new WstarR" << '\t' << i << std::endl;
								W(i, 0) = dshockR;
								W(i, 1) = ustar;
								W(i, 2) = pstar;

							}
						}
						//Right Rarefraction fan
						else {
							//Right of fastest rarefraction wave (unaffected by rarefraction)
							if (S > SHR){
								std::cout << "rarefraction WstarR" << '\t' << i << std::endl;
								W.row(i) = WstarR;
							}
							//Out of the rarefraction fan, within the right star state
							else if (S < STR){
								std::cout << "rarefraction new WstarR" << '\t' << i << std::endl;
								W(i, 0) = drareR;
								W(i, 1) = ustar;
								W(i, 2) = pstar;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								std::cout << "rarefraction WstarR" << '\t' << i << std::endl;
								double dRfan = dR*pow(C(2, 5) - (C(2, 6)/soundspeedR)*(uR - S), C(2, 4));
								double uRfan = C(2, 5)*(-soundspeedR + C(2, 7)*uR + S);
								double pRfan = pR*pow(C(2, 5) - (C(2, 6)/soundspeedR)*(uR - S), C(2, 3));
								vector WRfan(dRfan, uRfan, pRfan);
								W.row(i) = WRfan;
							}
						}
					}	
				}



/*----------------------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------------------*/

				if (S2 < ustar_2){ //S2 < SL_2){
					std::cout << "case 1 " << std::endl;
					double S = (xPos - Test.x0)/Test.tstop;
					if (S1 <= ustar_1){
						//Left Shock wave
						if (pstar_1 > pL_1){
							//Left of Shock Wave (unaffected by shock)
							if (S < SL_1){
								W.row(i) = WL;
								std::cout << "Shock WL" << '\t' << i << std::endl;
							}
							//Right of Shock Wave (shocked material, star state)
							else {
								W(i, 0) = dshockL_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;
								std::cout << "Shock WL star" << '\t' << i << std::endl;
							}
						}
						//Left Rarefraction fan
						else {
							//Left of fastest rarefraction wave (unaffected by rarefraction)
							if (S < SHL_1){
								W.row(i) = WL;
								std::cout << "Rarefraction WL" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the left star state
							else if (S > STL_1){
								W(i, 0) = drareL_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;
								std::cout << "Rarefraction WL star" << '\t' << i << std::endl;
							}
							//Within the rarefraction fan
							else {
								double dLfan = dL_1*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL_1 - S), C(0, 4));
								double uLfan = C(0, 5)*(C(0, 0) + C(0, 7)*uL_1 + S);
								double pLfan = pL_1*pow(C(0, 5) + (C(0, 6)/C(0, 0))*(uL_1 - S), C(0, 3));
								vector WLfan(dLfan, uLfan, pLfan);
								W.row(i) = WLfan;
								std::cout << "Rarefraction WL fan" << '\t' << i << std::endl;
							}
						}
					}

					else {
						//Right Shock wave
						if (pstar_1 > pR_1){
							//Right of Shock Wave (unaffected by shock)
							if (S > SR_1){
								std::cout << "Shock WM" << '\t' << i << std::endl;
								W.row(i) = WM1;

							}
							//Left of Shock Wave (shocked material, star state)
							else {
								std::cout << "Shock WM Star" << '\t' << i << std::endl;
								W(i, 0) = dshockR_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;

							}
						}
						//Right Rarefraction fan
						else {
							//Right of fastest rarefraction wave (unaffected by rarefraction)
							if (S > SHR_1){
								std::cout << "rarefraction WM" << '\t' << i << std::endl;
								W.row(i) = WM1;
							}
							//Out of the rarefraction fan, within the right star state
							else if (S < STR_1){
								std::cout << "rarefraction WM Star" << '\t' << i << std::endl;
								W(i, 0) = drareR_1;
								W(i, 1) = ustar_1;
								W(i, 2) = pstar_1;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								std::cout << "rarefraction WM fan" << '\t' << i << std::endl;
								double dRfan = dR_1*pow(C(2, 5) - (C(2, 6)/C(2, 0))*(uR_1 - S), C(2, 4));
								double uRfan = C(2, 5)*(-C(2, 0) + C(2, 7)*uR_1 + S);
								double pRfan = pR_1*pow(C(2, 5) - (C(2, 6)/C(2, 0))*(uR_1 - S), C(2, 3));
								vector WRfan(dRfan, uRfan, pRfan);
								W.row(i) = WRfan;
							}
						}
					}
				}

				else if (S2 < ustar_2 && S2 > SL_2){
					std::cout << "case 2" << std::endl;
					double S = (xPos - Test.x1)/Test.tstop;
					if (S1 < ustar){
						//Left Shock wave
						if (pstar > pL){
							//Left of Shock Wave (unaffected by shock)
							if (S < SL){
								W.row(i) = WstarL;
								std::cout << "Shock WstarL" << '\t' << i << std::endl;
							}
							//Right of Shock Wave (shocked material, star state)
							else {
								W(i, 0) = dshockL;
								W(i, 1) = ustar;
								W(i, 2) = pstar;
								std::cout << "Shock new WstarL" << '\t' << i << std::endl;
							}
						}
						//Left Rarefraction fan
						else {
							//Left of fastest rarefraction wave (unaffected by rarefraction)
							if (S < SHL){
								W.row(i) = WstarL;
								std::cout << "Rarefraction WstarL" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the left star state
							else if (S > STL){
								W(i, 0) = drareL;
								W(i, 1) = ustar;
								W(i, 2) = pstar;
								std::cout << "Rarefraction new WstarL" << '\t' << i << std::endl;
							}
							//Within the rarefraction fan
							else {
								double dLfan = dL*pow(C(2, 5) + (C(2, 6)/soundspeedL)*(uL - S), C(2, 4));
								double uLfan = C(2, 5)*(soundspeedL + C(2, 7)*uL + S);
								double pLfan = pL*pow(C(2, 5) + (C(2, 6)/soundspeedL)*(uL - S), C(2, 3));
								vector WLfan(dLfan, uLfan, pLfan);
								W.row(i) = WLfan;
								std::cout << "Rarefraction new WfanL" << '\t' << i << std::endl;
							}
						}
					}

					else {
						//Right Shock wave
						if (pstar > pR){
							//Right of Shock Wave (unaffected by shock)
							if (S > SR){
								std::cout << "Shock WstarR" << '\t' << i << std::endl;
								W.row(i) = WstarR;

							}
							//Left of Shock Wave (shocked material, star state)
							else {
								std::cout << "Shock new WstarR" << '\t' << i << std::endl;
								W(i, 0) = dshockR;
								W(i, 1) = ustar;
								W(i, 2) = pstar;

							}
						}
						//Right Rarefraction fan
						else {
							//Right of fastest rarefraction wave (unaffected by rarefraction)
							if (S > SHR){
								std::cout << "rarefraction WstarR" << '\t' << i << std::endl;
								W.row(i) = WstarR;
							}
							//Out of the rarefraction fan, within the right star state
							else if (S < STR){
								std::cout << "rarefraction new WstarR" << '\t' << i << std::endl;
								W(i, 0) = drareR;
								W(i, 1) = ustar;
								W(i, 2) = pstar;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								std::cout << "rarefraction WstarR" << '\t' << i << std::endl;
								double dRfan = dR*pow(C(2, 5) - (C(2, 6)/soundspeedR)*(uR - S), C(2, 4));
								double uRfan = C(2, 5)*(-soundspeedR + C(2, 7)*uR + S);
								double pRfan = pR*pow(C(2, 5) - (C(2, 6)/soundspeedR)*(uR - S), C(2, 3));
								vector WRfan(dRfan, uRfan, pRfan);
								W.row(i) = WRfan;
							}
						}
					}
				}

				else {
					double S = (xPos - Test.x1)/Test.tstop;

					//Left side of Contact wave
					if (S <= ustar_2){
						//Left Shock wave
						if (pstar_2 > pL_2){
							//Left of Shock Wave (unaffected by shock)
							if (S < SL_2){
								W.row(i) = WM1;
								std::cout << "right Shock WM" << '\t' << i << std::endl;
							}
							//Right of Shock Wave (shocked material, star state)
							else {
								W(i, 0) = dshockL_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "right Shock WM Star" << '\t' << i << std::endl;
							}
						}
						//Left Rarefraction fan
						else {
							//Left of fastest rarefraction wave (unaffected by rarefraction)
							if (S < SHL_2){
								W.row(i) = WM1;
								std::cout << "right rarefraction WM" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the left star state
							else if (S > STL_2){
								W(i, 0) = drareL_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "right rarefraction WM Star" << '\t' << i << std::endl;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								double dLfan = dL_2*pow(C(2, 5) + (C(2, 6)/C(2, 0))*(uL_2 - S), C(2, 4));
								double uLfan = C(2, 5)*(C(2, 0) + C(2, 7)*uL_2 + S);
								double pLfan = pL_2*pow(C(2, 5) + (C(2, 6)/C(2, 0))*(uL_2 - S), C(2, 3));
								vector WLfan(dLfan, uLfan, pLfan);
								W.row(i) = WLfan;
								std::cout << "right rarefraction WM fan" << '\t' << i << std::endl;
							}
						}
					}
					//right side of contact wave
					else {
						double S = (xPos - Test.x1)/Test.tstop;
						//Right Shock wave
						if (pstar_2 > pR_2){
							//Right of Shock Wave (unaffected by shock)
							if (S > SR_2){
								//std::cout << "Shock WR" << '\t' << i << std::endl;
								W.row(i) = WR;
								std::cout << "Shock WR" << '\t' << i << std::endl;
							}
							//Left of Shock Wave (shocked material, star state)
							else {
								//std::cout << "Shock Star" << '\t' << i << std::endl;
								W(i, 0) = dshockR_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "Shock WR Star" << '\t' << i << std::endl;
							}
						}
						//Right Rarefraction fan
						else {
							//Right of fastest rarefraction wave (unaffected by rarefraction)
							if (S > SHR_2){
								//std::cout << "rarefraction WR" << '\t' << i << std::endl;
								W.row(i) = WR;
								std::cout << "rarefraction WR" << '\t' << i << std::endl;
							}
							//Out of the rarefraction fan, within the right star state
							else if (S < STR_2){
								//std::cout << "rarefraction Star" << '\t' << i << std::endl;
								W(i, 0) = drareR_2;
								W(i, 1) = ustar_2;
								W(i, 2) = pstar_2;
								std::cout << "rarefraction WR Star" << '\t' << i << std::endl;
								//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
							}
							//Within the rarefraction fan
							else {
								//std::cout << "rarefraction fan" << '\t' << i << std::endl;
								double dRfan = dR_2*pow(C(1, 5) - (C(1, 6)/C(1, 0))*(uR_2 - S), C(1, 4));
								double uRfan = C(1, 5)*(-C(1, 0) + C(1, 7)*uR_2 + S);
								double pRfan = pR_2*pow(C(1, 5) - (C(1, 6)/C(1, 0))*(uR_2 - S), C(1, 3));
								vector WRfan(dRfan, uRfan, pRfan);
								W.row(i) = WRfan;
								std::cout << "rarefraction WR fan" << '\t' << i << std::endl;
							}
						}
					}	
				}
			}