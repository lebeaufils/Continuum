#include "Solvers.h"

/*--------------------------------------------------------------------------------
 * MUSCL
 --------------------------------------------------------------------------------*/
MUSCL::MUSCL(gfmTests Test)
	:CFL(0), N(Test.N), count(0), dt(0), dx(Test.L/Test.N), X(Test.N+2, 1), U(Test.N+4, 3), F(Test.N+2, 3), Smax(0), ULi(Test.N+4, 3), URi(Test.N+4, 3){
}

MUSCL::MUSCL(double c, eulerTests Test)
	:CFL(c), N(Test.N), count(0), dt(0), dx(Test.L/Test.N), X(Test.N+2, 1), U(Test.N+4, 3), F(Test.N+2, 3), Smax(0), ULi(Test.N+4, 3), URi(Test.N+4, 3){
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

	for (int i=0; i<N; i++){
		X(i+1) = (i+0.5)*dx;
		if (X(i+1)  < Test.x0){
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
	flux(1) = U(1)*(U(1)/U(0)) + (IG->y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	flux(2) = (U(1)/U(0))*(U(2) + (IG->y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)));
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
	/*	
		double Ppvrs, p0;
		double Pstar;
		double ql, qr;
		//------------------------------------------------------
		double rhoavg = 0.5*(dl + dr);
		double aavg = 0.5*(al + ar);
		Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg; //initial pressure estimate
		//------------------------------------------------------
		double Pmax = fmax(Pl, Pr);
		double Pmin = fmin(Pl, Pr);
		double Quser = 2.0; //Toro sect 9.5.1
		double z = (IG->y - 1)/(2*IG->y);

		//if the pressure difference is small and the guess is within the initial pressure range
		if (Pmax/Pmin < Quser && Ppvrs > Pmin && Ppvrs < Pmax){
			//Select the linearised pressure guess
			p0 = Ppvrs;
		}
		//If the guess is less than the initial pressure, two rarefraction wves have been formed
		else if (Ppvrs < Pmin){
			double pTR = pow((al + ar - 0.5*(IG->y - 1)*(ur - ul))/((al/pow(Pl, z)) + (ar/pow(Pr, z))), 1./z);
			p0 = pTR;
		}
		//If the pressure difference is large or if the guess is larger than Pmax
		else {
			//Select the two shock initial guess using pPV as estimate
			double AL = (2./(IG->y + 1))/dl;	double AR = (2./(IG->y + 1))/dr;
			double BL = Pl*((IG->y - 1)/(IG->y + 1));	double BR = Pr*((IG->y - 1)/(IG->y + 1));
			double QL = sqrt(AL/(Ppvrs + BL));
			double QR = sqrt(AR/(Ppvrs + BR));
			double pTS = (QL*Pl + QR*Pr - (ur - ul))/(QL + QR);
			p0 = pTS;
		}

		Pstar = fmax(0.0, p0);
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

		double SRpressure = ur + ar*qr;
		double SLpressure = ul - al*ql;
	*/	
		

		/*---------------------------------------
		 * Direct Wave Speed estimates
		 ---------------------------------------*/
		//Roe Riemann Solver
		double ubar = (sqrt(dl)*ul + sqrt(dr)*ur)/(sqrt(dl) + sqrt(dr));
		//double Hl = (hllcUL(2) + Pl)/dl; //enthalpy
		//double Hr = (hllcUR(2) + Pr)/dr;
		//double Hbar = (sqrt(dl)*Hl + sqrt(dr)*Hr)/(sqrt(dl) + sqrt(dr));
		//double abar = pow((IG->y - 1)*(Hbar - 0.5*pow(ubar, 2)), 0.5);

		//HLLE
		double n2 = 0.5*sqrt(dl)*sqrt(dr)/pow((sqrt(dl) + sqrt(dr)), 2);
		double d2 = (sqrt(dl)*pow(al, 2) + sqrt(dr)*pow(ar, 2))/(sqrt(dl) + sqrt(dr)) + n2*pow((ur - ul), 2);

		//Davies
		//double Splus = fmax(abs(ul) + al, abs(ur) + ar);

		//double SLdavies = -Splus; 			double SRdavies = -Splus;
		double SLhlle = ubar - sqrt(d2); 	double SRhlle = ubar + sqrt(d2);
		//double SLroe = ubar - abar; 		double SRroe = ubar + abar;
		SL = SLhlle; SR = SRhlle;
		
		if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

		Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));
		//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;
		//initialize FL and FR for each timestep
		vector FL(mvl, mvl*ul + Pl, ul*(El + Pl));
		vector FR(mvr, mvr*ur + Pr, ur*(Er + Pr));
		

		/*---------------------------------------
		 * davies wave speed estimate
		 ---------------------------------------*/

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
		//std::cout << FL << '\t' << FR << std::endl;
}

void MUSCL::conservative_update_formula(int i){
	U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
}

void MUSCL::conservative_update_formula(double newdt, double newdx, int i){
	U.row(i) = U.row(i) - (newdt/newdx)*(F.row(i) - F.row(i-1));
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
		std::cout << Smax << std::endl;
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

	for (int i=2; i<N+2; i++){

		double u = U(i, 1)/U(i, 0);
		double P = IG->Pressure(U, i);
		double e = IG->internalE(U, i);

		outfile << X(i-1) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	outfile.close();
	std::cout << "done: MUSCL" << std::endl;
}


/*--------------------------------------------------------------------------------
 * Exact Solver
 --------------------------------------------------------------------------------*/
//Note, W contains primitive variables (density, velocity, pressure)

EXACT::EXACT(int N, double dx, double x0, double y, vector WL, vector WR)
	: dx(dx), x0(x0), W(N, 3), WL(WL), WR(WR), TOL(1e-6), y(y), cL(0), cR(0), 
	CONST1(0),CONST2(0), CONST3(0), CONST4(0), CONST5(0), CONST6(0), CONST7(0), CONST8(0) {}

void EXACT::initial_conditions(){
	cL = sqrt(y*WL(2)/WL(0));
	cR = sqrt(y*WR(2)/WR(0));

	CONST1 = (y-1)/(2*y); // y-1 / 2y
	CONST2 = (y+1)/(2*y); // y+1 / 2y
	CONST3 = (2*y)/(y-1);
	CONST4 = 2./(y-1);
	CONST5 = 2./(y+1);
	CONST6 = (y-1)/(y+1);
	CONST7 = (y-1)/2.;
	CONST8 = y-1;
}

void EXACT::check_pressure_pos_condition(){
//(Δu)crit ≡ 2aL/γ−1 + 2aR/γ−1 ≤ uR −uL , (4.82) Toro pg 127
	double du_crit = CONST4*(cL + cR);
	double du = WR(1) - WL(1);

	//If the pressure positivity condition is not satisfied, vacuum is generated
	// and the solver fails. Exit the program if this condition is not satisfied.
	if (du_crit <= du){ //This conditions ensures S*L <= S*R
		throw "Pressure positivity condition violated";
	}
}

double EXACT::fk(double P, vector Wk){
	//data-dependent constants
	double Ak = CONST5/Wk(0);
	double Bk = Wk(2)*CONST6;
	double Qk = sqrt(Ak/(P + Bk)); //mass flux
	double ck = sqrt(y*Wk(2)/Wk(0));
	double flux;

	if (P > Wk(2)){ //Shock
		flux = (P - Wk(2))*Qk;
	}

	else { //Rarefraction
		flux = CONST4*ck*(pow(P/Wk(2), CONST1) - 1);
	}

	return flux;
}

double EXACT::f(double P){
	double du = (WR(1) - WL(1));
	return fk(P, WL) + fk(P, WR) + du;   
}

double EXACT::fkprime(double P, vector Wk){
	//data-dependent constants
	double Ak = CONST5/Wk(0);
	double Bk = Wk(2)*CONST6;
	double Qk = sqrt(Ak/(P + Bk)); //mass flux
	double ck = sqrt(y*Wk(2)/Wk(0));
	double fluxprime;

	if (P > Wk(2)){ //Shock
		fluxprime = Qk*(1 - (P - Wk(2))/(2.*(P + Bk))); 
	}

	else { //Rarefraction
		fluxprime = (1./(Wk(0)*ck))*pow(P/Wk(2), -CONST2);
	}

	return fluxprime;
}

double EXACT::fprime(double P){
	return fkprime(P, WR) + fkprime(P, WL);  
}

double EXACT::newton_raphson(double Pk){
	double Pk_1 = Pk - f(Pk)/fprime(Pk);
	return Pk_1;
}

double EXACT::relative_pressure_change(double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
	double CHA = 2*abs((Pk_1 - Pk)/(Pk_1 + Pk));
	return CHA;
}

double EXACT::compute_star_pressure(){
	//check positivity condition
	
	try {
		check_pressure_pos_condition();
	} 
	catch(const char* c){
		std::cout << c << std::endl;
		std::cout << "vacuum generated, terminating program" << std::endl;
	}
	
	//An approximation for p, p0 is required for the initial guess.
	//A poor choice of p0 results in the need for large number of iterations to achieve convergence
	double pPV, p0;
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
		double pTR = pow((cL + cR - 0.5*CONST8*(WR(1) - WL(1)))/((cL/pow(WL(2), CONST1)) + (cR/pow(WR(2), CONST1))), CONST3);
		p0 = pTR;
	}
	//If the pressure difference is large or if the guess is larger than Pmax
	else {
		//Select the two shock initial guess using pPV as estimate
		double AL = CONST5/WL(0);	double AR = CONST5/WR(0);
		double BL = WL(2)*CONST6;	double BR = WR(2)*CONST6;
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
}

double EXACT::compute_star_velocity(double pstar){
	double ustar = 0.5*(WL(1) + WR(1)) + 0.5*(fk(pstar, WR) - fk(pstar, WL));
	return ustar;
}

double EXACT::compute_shock_density(vector Wk, double pstar){
	//turn this into a stored variable so the iteration does not have to be called multiple times
	
	//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
	double Pratio = pstar/Wk(2);
	double dshock_k = Wk(0)*((CONST6 + Pratio)/(CONST6*Pratio + 1));
	return dshock_k;
}

double EXACT::compute_rarefraction_density(vector Wk, double pstar){
	double drare_k = Wk(0)*pow(pstar/Wk(2), 1./y);
	return drare_k;
}

void EXACT::sampling(double t){
	//Solution of wavestructure within the domain of interest
	double pstar = compute_star_pressure();
	double ustar = compute_star_velocity(pstar);
	double dshockL = compute_shock_density(WL, pstar);
	double dshockR = compute_shock_density(WR, pstar);
	double drareL = compute_rarefraction_density(WL, pstar);
	double drareR = compute_rarefraction_density(WR, pstar);

	//primitive variables of left and right states
	double dL = WL(0); double dR = WR(0);
	double uL = WL(1); double uR = WR(1);
	double pL = WL(2); double pR = WR(2);

	//Shock Speeds
	double SL = uL - cL*sqrt((CONST2*(pstar/pL) + CONST1)); 
	double SR = uR + cR*sqrt((CONST2*(pstar/pR) + CONST1));

	//Rarefraction wave speeds
	double cstarL = cL*pow(pstar/pL, CONST1);
	double cstarR = cR*pow(pstar/pR, CONST1);

	double SHL = uL - cL; double SHR = uR + cR; //Head of fan
	double STL = ustar - cstarL; double STR = ustar + cstarR; //tail of fan

	std::cout << drareL << '\t' << drareR << std::endl;
	std::cout << SHL << '\t' << SHR << std::endl;
	std::cout << STL << '\t' << STR<< std::endl;

	//Sampling is based on the position of wave with respect to time in x-t space,
	// characterised by its "speed" S = x/t.
	//When the solution at a specified time t is required the solution profiles are only a function of space x. Toro 137
	for (int i=0; i<N; i++){
		double xPos = (static_cast<double>(i) + 0.5)*dx;
		double S = (xPos - x0)/t;

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
					double dLfan = dL*pow(CONST5 + (CONST6/cL)*(uL - S), CONST4);
					double uLfan = CONST5*(cL + CONST7*uL + S);
					double pLfan = pL*pow(CONST5 + (CONST6/cL)*(uL - S), CONST3);
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
					W.row(i) = WR;
				}
				//Left of Shock Wave (shocked material, star state)
				else {
					W(i, 0) = dshockR;
					W(i, 1) = ustar;
					W(i, 2) = pstar;
				}
			}
			//Right Rarefraction fan
			else {
				//Right of fastest rarefraction wave (unaffected by rarefraction)
				if (S > SHR){
					W.row(i) = WR;
				}
				//Out of the rarefraction fan, within the right star state
				else if (S < STR){
					W(i, 0) = drareR;
					W(i, 1) = ustar;
					W(i, 2) = pstar;
					//std::cout << i << '\t' << S << '\t' << W.row(i) << std::endl;
				}
				//Within the rarefraction fan
				else {
					double dRfan = dR*pow(CONST5 - (CONST6/cR)*(uR - S), CONST4);
					double uRfan = CONST5*(-cR + CONST7*uR + S);
					double pRfan = pR*pow(CONST5 - (CONST6/cR)*(uR - S), CONST3);
					vector WRfan(dRfan, uRfan, pRfan);
					W.row(i) = WRfan;
				}
			}
		}
	}
}

void EXACT::output(){
	std::ofstream outfile;
	outfile.open("dataexact.txt");

	for (int i=0; i<N; i++){

		double xPos = (static_cast<double>(i) + 0.5)*dx;
		double e = W(i, 2)/(W(i, 0)*(y-1));

		outfile << xPos << '\t' << W(i, 0) << '\t' << W(i, 1)
				<< '\t' << W(i, 2) << '\t' << e << std::endl;
	}
	outfile.close();
	std::cout << "done: exact" << std::endl;
}

/*
double EXACT::fk(double P, vector Wk, EOS* IG){
	//Data-dependent constants
	double Ak = 2./((IG->y + 1)*Wk(0));
	double Bk = Wk(2)*((IG->y - 1)/(IG->y + 1));
	double ck = pow(IG->y*Wk(2)/Wk(0), 0.5); //soundspeed
	double flux;

	if (P > Wk(2)){ //Shock
		flux = (P - Wk(2))*pow(Ak/(P + Bk), 0.5);
	}

	else { //Rarefraction
		flux = ((2.*ck)/(IG->y - 1))*(pow(P/Wk(2), (IG->y - 1)/(2.*IG->y)) - 1);
	}

	return flux;
}

double EXACT::f(double P, vector WL, vector WR, EOS* IG){
	double du = (WR(1) - WL(1));
	return fk(P, WL, IG) + fk(P, WR, IG) + du;                                                                                                             
}

double EXACT::fkprime(double P, vector Wk, EOS* IG){ //first derivative of f
	//Data-dependent constants
	double Ak = 2./((IG->y + 1)*Wk(0));
	double Bk = Wk(2)*((IG->y - 1)/(IG->y + 1));
	double ck = pow(IG->y*Wk(2)/Wk(0), 0.5); //soundspeed
	double fluxprime;

	if (P > Wk(2)){ //Shock
		fluxprime = pow(Ak/(P + Bk), 0.5) * (1 - (P - Wk(2))/(2.*(P + Bk))); 
	}

	else { //Rarefraction
		fluxprime = (1./(Wk(0)*ck))*pow(P/Wk(2), -(IG->y + 1)/(2*IG->y));
	}

	return fluxprime;
}

double EXACT::fprime(double P, vector WL, vector WR, EOS* IG){
	return fkprime(P, WR, IG) + fkprime(P, WL, IG);
}

double EXACT::newton_raphson(double Pk, vector WL, vector WR, EOS* IG){
	double Pk_1 = Pk - f(Pk, WL, WR, IG)/fprime(Pk, WL, WR, IG);
	//double du = (WR(1) - WL(1));
	//double Pk_1 = Pk - (fk(Pk, WR, IG) + fk(Pk, WL, IG) + du)/(fkprime(Pk, WR, IG) + fkprime(Pk, WL, IG));
	return Pk_1;
}

double EXACT::relative_pressure_change(double Pk_1, double Pk){ //where Pk_1 is the k+1th iterate
	double CHA = abs(Pk_1 - Pk)/(0.5*(Pk_1 + Pk));
	return CHA;
}

double EXACT::compute_star_pressure(vector WL, vector WR, EOS* IG){
	//An approximation for p, p0 is required for the initial guess.
	//A poor choice of p0 results in the need for large number of iterations to achieve convergence
	
	double cL = pow(IG->y*WL(2)/WL(0), 0.5); //soundspeed in left state
	double cR = pow(IG->y*WR(2)/WR(0), 0.5); //soundspeed in right state

	double exponent = (IG->y - 1)/(2*IG->y);
	double _exponent = (2*IG->y)/(IG->y - 1);
	
	//Two-Rarefraction approximation
	double pTR = pow((cL + cR - 0.5*(IG->y-1)*(WR(1) - WL(1)))/((cL/pow(WL(2), exponent)) + (cR/pow(WR(2), exponent))), _exponent);
		//Pressure under the assumption that the two non-linear waves are rarefraction waves

	//double p0 = 0.5*(WL(2) + WR(2));
	//to make logic gates to select the 4 different P guesses...

	double Pk; Pk = pTR; //First guess
	double Pk_1;
	double CHA = 0;

	int count = 0;
	do{
		Pk_1 = newton_raphson(Pk, WL, WR, IG);
		CHA = relative_pressure_change(Pk_1, Pk);
		Pk = Pk_1; //Set the iterate as the new guess
		count += 1;
	}while(CHA > TOL);

	std::cout << count <<std::endl;
	return Pk;
}

double EXACT::compute_star_velocity(vector WL, vector WR, EOS* IG, double pstar){
	ustar = 0.5*(WL(1) + WR(1)) + 0.5*(fk(pstar, WR, IG) - fk(pstar, WL, IG));
	return ustar;
}

double EXACT::compute_star_density_k(vector Wk, EOS* IG, double pstar){
	//turn this into a stored variable so the iteration does not have to be called multiple times
	
	//From the Hugoniot jump conditions, see TORO 3.1.3 (substituting the expressio n for internal energy)
	Yratio = (IG.y-1)/(IG.y+1);
	Pratio = pstar/Wk(2);
	dstar_k = Wk(0)*((Yratio + Pratio)/(Yratio*Pratio + 1));
	return dstar_k;
}

double EXACT::compute_mass_flux_k(vector Wk, EOS* IG, double pstar){
	double Ak = 2./((IG->y + 1)*Wk(0));
	double Bk = Wk(2)*((IG->y - 1)/(IG->y + 1));
	Qk = pow((pstar + Bk)/Ak, 0.5);
	return Qk;
}

double EXACT::compute_left_shock_speed(vector WL, EOS* IG, double pstar){
	QL = compute_mass_flux_k(WL, IG, Pstar);
	return WL(1) - QL/WL(0);
}

double EXACT::compute_right_shock_speed(vector WR, EOS* IG, double pstar){
	QR = compute_mass_flux_k(WR, IG, Pstar);
	return WR(1) + QR/WR(0);
}
*/







