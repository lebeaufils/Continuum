#include "MUSCL.h"


MUSCL::MUSCL(double c, eulerTests Test)
	:CFL(c), X(Test.N+4, 1), U(Test.N+4, 3), F(Test.N+2, 3), count(0), N(Test.N), dx(Test.L/Test.N), dt(0){
}

void MUSCL::boundary_conditions(){
	U.row(1) = U.row(2); //muscl requires extra ghost cells on either boundary.
	U.row(0) = U.row(1);
	U.row(N+2) = U.row(N+1);
	U.row(N+3) = U.row(N+2);
}

void MUSCL::initial_conditions(IdealGas IG, eulerTests Test){
	vector eulerL;
	vector eulerR;

	//rho
	eulerL(0) = Test.initialL(0);
	eulerR(0) = Test.initialR(0);

	//rhou
	eulerL(1) = Test.initialL(0)*Test.initialL(1);
	eulerR(1) = Test.initialR(0)*Test.initialR(1);

	//E
	eulerL(2) = Test.initialL(2)/(IG.y-1) + 0.5*Test.initialL(0)*Test.initialL(1)*Test.initialL(1);
	eulerR(2) = Test.initialR(2)/(IG.y-1) + 0.5*Test.initialR(0)*Test.initialR(1)*Test.initialR(1);

	for (int i=0; i<N+2; i++){
		X(i+2) = i*dx;
		if (X(i+2)  < Test.x0){
			U.row(i+2) = eulerL;
		}
		else U.row(i+2) = eulerR;
	}

	boundary_conditions();

}

void MUSCL::initial_conditions(JWL MG, eulerTests Test){
	vector eulerL;
	vector eulerR;

	//rho
	eulerL(0) = Test.initialL(0);
	eulerR(0) = Test.initialR(0);

	//rhou
	eulerL(1) = Test.initialL(0)*Test.initialL(1);
	eulerR(1) = Test.initialR(0)*Test.initialR(1);

	//E
	eulerL(2) = Test.initialL(0)*(0.5*pow(Test.initialL(1), 2.0) + Test.initialL(2)/(Test.initialL(0)*MG.Gruneisen));
	eulerR(2) = Test.initialR(0)*(0.5*pow(Test.initialR(1), 2.0) + Test.initialR(2)/(Test.initialR(0)*MG.Gruneisen));

	for (int i=0; i<N+2; i++){
		X(i+2) = i*dx;
		if (X(i+2)  < Test.x0){
			U.row(i+2) = eulerL;
		}
		else U.row(i+2) = eulerR;
	}

	boundary_conditions();
}

vector MUSCL::f(vector U, IdealGas IG){
	vector flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + (IG.y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	flux(2) = (U(1)/U(0))*(U(2) + (IG.y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)));
	return flux;
}

vector MUSCL::f(vector U, JWL MG){
	vector flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + MG.PressureScalar(U);
	flux(2) = (U(1)/U(0))*(U(2) + MG.PressureScalar(U));
	return flux;
}

vector MUSCL::superBee(int i){
	//Calculating slope
	vector diMinus = U.row(i) - U.row(i-1);
	vector diPlus = U.row(i+1) - U.row(i);
	vector di = 0.5*diMinus + 0.5*diPlus;

	//calculating limiter
	vector epsiloni(1, 1, 1);
	vector epsilonR(0, 0, 0);
	vector ri(0, 0, 0);

	for (int j=0; j<3; j++){
		ri(j) = diMinus(j)/diPlus(j);
	}

	for (int j=0; j<3; j++){
		if (ri(j) <= 0){
			epsiloni(j) = 0;
		}

		else if (0 < ri(j) && ri(j) <= 0.5){
			epsiloni(j) = 2*ri(j);
		}

		else if (0.5 < ri(j) && ri(j) <= 1){
			epsiloni(j) = 1.0;
		}

		else if (ri(j) >= 1){
			//epsilonR(j) = 2*(2/(1-ctmp(j)))/(1+ri(j));
			epsilonR(j) = 2/(1+ri(j));
			epsiloni(j) = fmin(fmin(ri(j), epsilonR(j)), 2.0);
		}
	}

	vector diBar(0, 0, 0);
	for (int j=0; j<3; j++){
		 diBar(j) = epsiloni(j)*di(j);
	}

	return diBar;
}

vector MUSCL::vanLeer(int i){

	vector diMinus = U.row(i) - U.row(i-1);
	vector diPlus = U.row(i+1) - U.row(i);
	vector di = 0.5*diMinus + 0.5*diPlus;

	/*-----------------------------------------------
	 * Slope limiter -- Van Leer
	 ----------------------------------------------*/

	vector epsiloni(1, 1, 1);
	vector epsilonR(0, 0, 0);
	vector ri(0, 0, 0);

	for (int j=0; j<3; j++){
		ri(j) = diMinus(j)/diPlus(j);
	}

	for (int j=0; j<3; j++){
		if (ri(j) <= 0){
			epsiloni(j) = 0;
		}
		//note: A more refined approach would be to adopt characteristic limiting p510
		else if (ri(j) > 0){
			epsilonR(j) = 2/(1+ri(j));
			epsiloni(j) = fmin(2*ri(j)/(1+ri(j)), epsilonR(j));
		}
	}

	vector diBar(0, 0, 0);
	for (int j=0; j<3; j++){
		 diBar(j) = epsiloni(j)*di(j);
	}

	return diBar;
	//make r = 0 to test first order
}

vector MUSCL::minBee(int i){


	vector diMinus = U.row(i) - U.row(i-1);
	vector diPlus = U.row(i+1) - U.row(i);
	vector di = 0.5*diMinus + 0.5*diPlus;

	/*-----------------------------------------------
	 * Slope limiter -- MinBee
	 ----------------------------------------------*/

	vector epsiloni(1, 1, 1);
	vector epsilonR(0, 0, 0);
	vector ri(0, 0, 0);

	for (int j=0; j<3; j++){
		ri(j) = diMinus(j)/diPlus(j);
	}

	for (int j=0; j<3; j++){
		if (ri(j) <= 0){
			epsiloni(j) = 0;
		}
		else if (0 < ri(j) &&  ri(j) <= 1){
			epsiloni(j) = ri(j);
		}
		else if (ri(j) > 1){
			epsilonR(j) = 2/(1+ri(j));
			epsiloni(j) = fmin(1.0, epsilonR(j));
		}
	}

	vector diBar(0, 0, 0);
	for (int j=0; j<3; j++){
		 diBar(j) = epsiloni(j)*di(j);
	}

	return diBar;
}

enum slopeLimiter {MinBee, VanLeer, SuperBee, Quit};
slopeLimiter getLimiter(){
	//static std::map<char, slopeLimiter> theMap = {{"M", MinBee}, {"V", VanLeer}, {"S", SuperBee}};

	static std::map<std::string, slopeLimiter> theMap;
		theMap["m"] = MinBee;
		theMap["v"] = VanLeer;
		theMap["s"] = SuperBee;
		theMap["q"] = Quit;

	std::string str;
	std::cout << "Slope limiter options" << std::endl
			<< "Min Bee (M)" << std::endl
			<< "Van Leer (V)" << std::endl
			<< "Super Bee (S)" << std::endl
			<< "enter (Q) to exit" << std::endl;
	do{
		std::cin >> str;
		bool find = theMap.count(str);

		if (find == 1){
			slopeLimiter a = theMap[str];
			switch(a){
			case MinBee:
				std::cout << "Min Bee slope-limiter" << std::endl;
				break;
			case VanLeer:
				std::cout << "Van Leer slope-limiter" << std::endl;
				break;
			case SuperBee:
				std::cout << "Super Bee slope-limiter" << std::endl;
				break;
			case Quit:
				exit(0);
			}
			break;
		}
		else{
			std::cout << "Invalid input, enter (M), (V) or (S). (Q) to exit." << std::endl;
		}

	}while(true);
	return theMap[str];
}

void MUSCL::solver(IdealGas IG, eulerTests Test){

	double al, ar;

	double ul, ur; //note u is the eigenvalue of the euler equation
	double dl, dr;
	double Pl, Pr;
	//Pressure-based wave speed estimate
	double Ppvrs;
	double Pstar;
	double rhoavg;
	double aavg;
	double ql, qr;

	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;

	double Smax = 0;
	vector tmp;

	//MUSCL
	matrix ULi; ULi.resize(N+4, 3);
	matrix URi; URi.resize(N+4, 3);
	//matrix Ftmp; Ftmp.resize(N+2, 3);
	//matrix Utmp; Utmp.resize(N+4, 3);
	vector Utmp;

	slopeLimiter a = getLimiter();

	double t = 0.0;
	do{

		for (int i=1; i<N+3; i++){ //U goes from 0 to N+3

			/*-----------------------------------------------
			 * Data Reconstruction
			 ----------------------------------------------*/
			switch(a){
			case MinBee:
				Utmp = U.row(i);
				ULi.row(i) = Utmp - 0.5*minBee(i);
				URi.row(i) = Utmp + 0.5*minBee(i);
				break;
			case VanLeer:
				Utmp = U.row(i);
				ULi.row(i) = Utmp - 0.5*vanLeer(i);
				URi.row(i) = Utmp + 0.5*vanLeer(i);
				break;
			case SuperBee:
				Utmp = U.row(i);
				ULi.row(i) = Utmp - 0.5*superBee(i);
				URi.row(i) = Utmp + 0.5*superBee(i);
				break;
			case Quit:
				exit(0);
			}

		}

			/*-----------------------------------------------
			 * Evolution by 1/2 time-step
			 ----------------------------------------------*/
		for (int i=1; i<N+2; i++){

			vector ULtmp = ULi.row(i);
			vector URtmp = URi.row(i);
			vector ULtmp1 = ULi.row(i+1);
			vector URtmp1 = URi.row(i+1);

			vector ULbar = ULtmp1 + 0.5*(dt/dx)*(f(ULtmp1, IG) - f(URtmp1, IG)); //UL(i+1)
			vector URbar = URtmp + 0.5*(dt/dx)*(f(ULtmp, IG) - f(URtmp, IG));


			/*-------------------------------------------------------
			 * Solution of the piecewise constant Riemann problem pg 180
			 -------------------------------------------------------*/
			/*-------------------------------------------------------
			 * HLLC solver
			 -------------------------------------------------------*/
			vector hllcUL = URbar;
			vector hllcUR = ULbar;

			//if (count == 1) std::cout << U.row(i) << std::endl;

			//conservative variables
				//density
				dl = hllcUL(0);
				dr = hllcUR(0);
				//momentum
				mvl = hllcUL(1);
				mvr = hllcUR(1);
				//energy
				El = hllcUL(2);
				Er = hllcUR(2);

				//Pressure
				Pl = (IG.y-1)*(hllcUL(2) - 0.5*hllcUL(0)*pow((hllcUL(1)/hllcUL(0)),2.0));
				Pr = (IG.y-1)*(hllcUR(2) - 0.5*hllcUR(0)*pow((hllcUR(1)/hllcUR(0)),2.0));

				//velocity
				ul = hllcUL(1)/hllcUL(0);
				ur = hllcUR(1)/hllcUR(0);

				//soundspeed
				al = sqrt(IG.y*(Pl/dl));
				ar = sqrt(IG.y*(Pr/dr));


				/*---------------------------------------
				 * pressure based wave speed estimate
				 ---------------------------------------*/

				rhoavg = 0.5*(dl + dr);
				aavg = 0.5*(al + ar);

				Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;

				Pstar = fmax(0.0, Ppvrs);
				if (Pstar <= Pl){
					ql = 1.0;
				}

				else if (Pstar > Pl){
					ql = sqrt(1 + ((IG.y+1)/(2*IG.y))*((Pstar/Pl) - 1));
				}

				if (Pstar <= Pr){
					qr = 1.0;
				}

				else if (Pstar > Pr){
					qr = sqrt(1 + ((IG.y+1)/(2*IG.y))*((Pstar/Pr) - 1));
				}

				SR = ur + ar*qr;
				SL = ul - al*ql;
				//finding Smax for the whole domain (for each timestep)
				if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

				Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));
				//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;
				//initialize FL and FR for each timestep
				vector FL(mvl, mvl*ul + Pl, ul*(El + Pl));
				vector FR(mvr, mvr*ur + Pr, ur*(Er + Pr));

				if (0 <= SL){
					F.row(i) = FL;
				}

				//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
				else if (SL<=0 && Sstar>=0){
					double tmpUstar = dl*((SL - ul)/(SL - Sstar));
					vector UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))));
					tmp = U.row(i);
					F.row(i) = FL + SL*(UstarL - tmp);
				}

				else if (Sstar<=0 && SR>=0){
					double tmpUstar = dr*((SR - ur)/(SR - Sstar));
					vector UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))));
					tmp = U.row(i+1);
					F.row(i) = FR + SR*(UstarR - tmp);
				}

				else if (0 >= SR){
					F.row(i) = FR;
				}
				//if (count == 0) std::cout << F.row(i) << std::endl;
		}
		//end of domain loop

		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > Test.tstop) dt = Test.tstop - t;
		t += dt;
		count += 1;

		if (count == 0) std::cout << dt << std::endl;
		//updating U
		for (int i=2; i<N+2; i++){
			//if (count == 0) std::cout << U.row(i) << '\t' << '\t';
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
			//if (count == 0) std::cout << U.row(i) << std::endl;
		}
		boundary_conditions();


	}while (t<Test.tstop);
	std::cout << count << std::endl;
}

void MUSCL::solver(JWL MG, eulerTests Test){

	double al, ar;

	double ul, ur; //note u is the eigenvalue of the euler equation
	double dl, dr;
	double Pl, Pr;

	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;
	double Splus, Smax = 0;

	vector tmp;

	//MUSCL
	matrix ULi; ULi.resize(N+4, 3);
	matrix URi; URi.resize(N+4, 3);

	vector Utmp;

	slopeLimiter a = getLimiter();

	double t = 0.0;
	do{

		for (int i=1; i<N+3; i++){ //U goes from 0 to N+3

			/*-----------------------------------------------
			 * Data Reconstruction
			 ----------------------------------------------*/
			switch(a){
			case MinBee:
				Utmp = U.row(i);
				ULi.row(i) = Utmp - 0.5*minBee(i);
				URi.row(i) = Utmp + 0.5*minBee(i);
				break;
			case VanLeer:
				Utmp = U.row(i);
				ULi.row(i) = Utmp - 0.5*vanLeer(i);
				URi.row(i) = Utmp + 0.5*vanLeer(i);
				break;
			case SuperBee:
				Utmp = U.row(i);
				ULi.row(i) = Utmp - 0.5*superBee(i);
				URi.row(i) = Utmp + 0.5*superBee(i);
				break;
			case Quit:
				exit(0);
			}

		}

			/*-----------------------------------------------
			 * Evolution by 1/2 time-step
			 ----------------------------------------------*/
		for (int i=1; i<N+2; i++){

			vector ULtmp = ULi.row(i);
			vector URtmp = URi.row(i);
			vector ULtmp1 = ULi.row(i+1);
			vector URtmp1 = URi.row(i+1);

			vector ULbar = ULtmp1 + 0.5*(dt/dx)*(f(ULtmp1, MG) - f(URtmp1, MG)); //UL(i+1)
			vector URbar = URtmp + 0.5*(dt/dx)*(f(ULtmp, MG) - f(URtmp, MG));


			/*-------------------------------------------------------
			 * Solution of the piecewise constant Riemann problem pg 180
			 -------------------------------------------------------*/
			/*-------------------------------------------------------
			 * HLLC solver
			 -------------------------------------------------------*/
			vector hllcUL = URbar;
			vector hllcUR = ULbar;

			//if (count == 1) std::cout << U.row(i) << std::endl;

			//conservative variables
				//density
				dl = hllcUL(0);
				dr = hllcUR(0);
				//momentum
				mvl = hllcUL(1);
				mvr = hllcUR(1);
				//energy
				El = hllcUL(2);
				Er = hllcUR(2);

				//Pressure
				Pl = MG.PressureScalar(hllcUL);
				Pr = MG.PressureScalar(hllcUR);

				//velocity
				ul = hllcUL(1)/hllcUL(0);
				ur = hllcUR(1)/hllcUR(0);

				//soundspeed
				al = MG.soundspeedScalar(hllcUL);
				ar = MG.soundspeedScalar(hllcUR);


				/*---------------------------------------
				 * Davies wave speed estimate
				 ---------------------------------------*/
				//S+ = max{| uL | +aL,| uR | +aR} .
				Splus = std::max(abs(ul) + al, abs(ur) + ar);
				SL = -Splus;
				SR = Splus;

				//finding Smax for the whole domain (for each timestep)
				if (Splus > Smax) Smax = Splus;

				Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));
				//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;
				//initialize FL and FR for each timestep
				vector FL(mvl, mvl*ul + Pl, ul*(El + Pl));
				vector FR(mvr, mvr*ur + Pr, ur*(Er + Pr));

				if (0 <= SL){
					F.row(i) = FL;
				}

				//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
				else if (SL<=0 && Sstar>=0){
					double tmpUstar = dl*((SL - ul)/(SL - Sstar));
					vector UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))));
					tmp = U.row(i);
					F.row(i) = FL + SL*(UstarL - tmp);
				}

				else if (Sstar<=0 && SR>=0){
					double tmpUstar = dr*((SR - ur)/(SR - Sstar));
					vector UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))));
					tmp = U.row(i+1);
					F.row(i) = FR + SR*(UstarR - tmp);
				}

				else if (0 >= SR){
					F.row(i) = FR;
				}
				//if (count == 0) std::cout << F.row(i) << std::endl;
		}
		//end of domain loop

		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > Test.tstop) dt = Test.tstop - t;
		t += dt;
		count += 1;

		if (count == 0) std::cout << dt << std::endl;
		//updating U
		for (int i=2; i<N+2; i++){
			//if (count == 0) std::cout << U.row(i) << '\t' << '\t';
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
			//if (count == 0) std::cout << U.row(i) << std::endl;
		}
		boundary_conditions();


	}while (t<Test.tstop);
	std::cout << count << std::endl;
}

void MUSCL::output(IdealGas IG){

	std::ofstream outfile;
	outfile.open("dataMUSCL.txt");

	for (int i=2; i<N+3; i++){

		double u = U(i, 1)/U(i, 0);
		double P = IG.Pressure(U, i);
		double e = IG.internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}

void MUSCL::output(JWL MG){

	std::ofstream outfile;
	outfile.open("dataMUSCL.txt");

	for (int i=2; i<N+3; i++){

		double u = U(i, 1)/U(i, 0);
		double P = MG.Pressure(U, i);
		double e = MG.internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}
