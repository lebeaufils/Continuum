/*
 * Solvers.C
 *
 *  Created on: 24 Jan 2019
 *      Author: forte
 */

#include "Solvers.h"

RPsolvers::RPsolvers(gfmTests Test, int nU, int nF)
	: CFL(0), N(Test.N), count(0), dt(0), dx(Test.L/Test.N), X(nU, 1), U(nU, 3), F(nF, 3), Smax(0){
}

RPsolvers::RPsolvers(double c, eulerTests Test, int nU, int nF)
	: CFL(c), N(Test.N), count(0), dt(0), dx(Test.L/Test.N), X(nU, 1), U(nU, 3), F(nF, 3), Smax(0){
}

void RPsolvers::conservative_update_formula(int i){
	U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
}

void RPsolvers::conservative_update_formula(double newdt, double newdx, int i){
	U.row(i) = U.row(i) - (newdt/newdx)*(F.row(i) - F.row(i-1));
}

/*--------------------------------------------------------------------------------
 * HLLC
 --------------------------------------------------------------------------------*/
HLLC::HLLC(gfmTests Test)
	: RPsolvers(Test, Test.N+2, Test.N+1){
}

HLLC::HLLC(double c, eulerTests Test)
	: RPsolvers(c, Test, Test.N+2, Test.N+1){
}

void HLLC::boundary_conditions(){
	U.row(0) = U.row(1);
	U.row(N+1) = U.row(N);
}

void HLLC::initial_conditions(EOS* IG, eulerTests Test){
	vector eulerL;
	vector eulerR;

	eulerL = IG->conservedVar(Test.initialL);
	eulerR = IG->conservedVar(Test.initialR);

	for (int i=0; i<N; i++){
		X(i+1) = i*dx;
		if (X(i+1)  < Test.x0){
			U.row(i+1) = eulerL;
		}
		else U.row(i+1) = eulerR;
	}

	boundary_conditions();
}

void HLLC::compute_fluxes(EOS* IG, int i){

	//Storing left and right states
	double al, ar; //sound-speed
	double ul, ur; //wave velocity
	double dl, dr; //density
	double Pl, Pr; //pressure

	//Pressure-based wave speed estimate
	//as in Toro's book chapter 10
	double Ppvrs;
	double Pstar;
	double rhoavg;
	double aavg;
	double ql, qr;

	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;

	vector tmp;

	//conservative variables
	//density
	dl = U(i, 0);
	dr = U(i+1, 0);
	//momentum
	mvl = U(i, 1);
	mvr = U(i+1, 1);
	//energy
	El = U(i, 2);
	Er = U(i+1, 2);

	//physical properties
	//Pressure
	Pl = IG->Pressure(U, i);
	Pr = IG->Pressure(U, i+1);

	//velocity
	ul = U(i, 1)/U(i, 0);
	ur = U(i+1, 1)/U(i+1, 0);

	//soundspeed
	al = IG->soundspeed(U, i);
	ar = IG->soundspeed(U, i+1);

	/*--------------------------------------------------------------
	 * pressure based wave speed estimate
	 *--------------------------------------------------------------*/

	rhoavg = 0.5*(dl + dr);
	aavg = 0.5*(al + ar);

	Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;
	Pstar = fmax(0.0, Ppvrs);

	if (Pstar <= Pl){
		ql = 1.0;
	}

	else{ //if (Pstar > Pl){
		ql = sqrt(1 + ((IG->y+1)/(2*IG->y))*((Pstar/Pl) - 1));
	}

	if (Pstar <= Pr){
		qr = 1.0;
	}

	else{ // if (Pstar > Pr){
		qr = sqrt(1 + ((IG->y+1)/(2*IG->y))*((Pstar/Pr) - 1));
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
	//Fluxes computed at each grid point
	//To be looped over domain grid
}


void HLLC::solver(EOS* IG, eulerTests Test){

	double t = 0.0;
	do{
		//compute fluxes at current timestep
		for (int i=0; i<N+1; i++){
			compute_fluxes(IG, i);
		}

		//set timestep following CFL conditions with max wavespeed Smax
		dt = CFL*(dx/Smax); 
		if (t + dt > Test.tstop) dt = Test.tstop - t; //

		//updating U
		for (int i=1; i<N+1; i++){
			//U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
			conservative_update_formula(i);
		}
		boundary_conditions();

		t += dt;
		count += 1;

	}while (t < Test.tstop);
	std::cout << count << std::endl;
}

void HLLC::output(EOS* IG){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	for (int i=1; i<N+1; i++){

		double u = U(i, 1)/U(i, 0);
		double P = IG->Pressure(U, i);
		double e = IG->internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done: HLLC" << std::endl;
}

/*--------------------------------------------------------------------------------
 * MUSCL
 --------------------------------------------------------------------------------*/
MUSCL::MUSCL(gfmTests Test)
	:RPsolvers(Test, Test.N+4, Test.N+2), ULi(Test.N+4, 3), URi(Test.N+4, 3){
}

MUSCL::MUSCL(double c, eulerTests Test)
	:RPsolvers(c, Test, Test.N+4, Test.N+2), ULi(Test.N+4, 3), URi(Test.N+4, 3){
}

void MUSCL::boundary_conditions(){
	U.row(1) = U.row(2); //muscl requires extra ghost cells on either boundary.
	U.row(0) = U.row(1);
	U.row(N+2) = U.row(N+1);
	U.row(N+3) = U.row(N+2);
}

void MUSCL::initial_conditions(EOS *IG, eulerTests Test){
	vector eulerL;
	vector eulerR;

	eulerL = IG->conservedVar(Test.initialL);
	eulerR = IG->conservedVar(Test.initialR);

/*
	//rho
	eulerL(0) = Test.initialL(0);
	eulerR(0) = Test.initialR(0);

	//rhou
	eulerL(1) = Test.initialL(0)*Test.initialL(1);
	eulerR(1) = Test.initialR(0)*Test.initialR(1);

	//E
	eulerL(2) = Test.initialL(2)/(IG->y-1) + 0.5*Test.initialL(0)*Test.initialL(1)*Test.initialL(1);
	eulerR(2) = Test.initialR(2)/(IG->y-1) + 0.5*Test.initialR(0)*Test.initialR(1)*Test.initialR(1);
*/

	for (int i=0; i<N+2; i++){
		X(i+2) = i*dx;
		if (X(i+2)  < Test.x0){
			U.row(i+2) = eulerL;
		}
		else U.row(i+2) = eulerR;
	}

	boundary_conditions();

}

/*
vector MUSCL::f(vector U, EOS IG){
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
*/
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

slopeLimiter MUSCL::getLimiter(){
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

//-----------------------------------------------------------------------------------

void MUSCL::data_reconstruction(slopeLimiter a){
	vector Utmp;

	//slopeLimiter a = getLimiter();

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
}

void MUSCL::compute_fluxes(EOS* IG, int i){

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

	vector tmp;


	/*-----------------------------------------------
	 * Evolution by 1/2 time-step
	 ----------------------------------------------*/

	vector ULtmp = ULi.row(i);
	vector URtmp = URi.row(i);
	vector ULtmp1 = ULi.row(i+1);
	vector URtmp1 = URi.row(i+1);

	//Post 1/2 time-step evolution left and right states
	vector ULbar = ULtmp1 + 0.5*(dt/dx)*(IG->f(ULtmp1) - IG->f(URtmp1)); //UL(i+1)
	vector URbar = URtmp + 0.5*(dt/dx)*(IG->f(ULtmp) - IG->f(URtmp));

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
		Pl = IG->PressureScalar(hllcUL);
		Pr = IG->PressureScalar(hllcUR);

		//velocity
		ul = hllcUL(1)/hllcUL(0);
		ur = hllcUR(1)/hllcUR(0);

		//soundspeed
		al = IG->soundspeedScalar(hllcUL);
		ar = IG->soundspeedScalar(hllcUR);


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

		else {
			ql = sqrt(1 + ((IG->y+1)/(2*IG->y))*((Pstar/Pl) - 1));
		}

		if (Pstar <= Pr){
			qr = 1.0;
		}

		else {
			qr = sqrt(1 + ((IG->y+1)/(2*IG->y))*((Pstar/Pr) - 1));
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

void MUSCL::solver(EOS* IG, eulerTests Test){
	
	slopeLimiter a = getLimiter();

	double t = 0.0;
	do{
		data_reconstruction(a);
		for (int i=1; i<N+2; i++){
			compute_fluxes(IG, i);
		}

		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > Test.tstop) dt = Test.tstop - t;
		t += dt;
		count += 1;

		if (count == 0) std::cout << dt << std::endl;
		//updating U
		for (int i=2; i<N+2; i++){
			conservative_update_formula(i);
		}
		boundary_conditions();

	}while (t < Test.tstop);
	std::cout << count << std::endl;
}

void MUSCL::output(EOS* IG){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	for (int i=2; i<N+3; i++){

		double u = U(i, 1)/U(i, 0);
		double P = IG->Pressure(U, i);
		double e = IG->internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done: MUSCL" << std::endl;
}









