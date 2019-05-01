#include "../headerFiles/Solvers.h"

/*--------------------------------------------------------------------------------
 * MUSCL
 --------------------------------------------------------------------------------*/

void MUSCL::boundary_conditions(Euler1D &var, Domain1D domain){
	var.U.row(1) = var.U.row(2); //muscl requires extra ghost cells on either boundary.
	var.U.row(0) = var.U.row(1);
	var.U.row(domain.N+2) = var.U.row(domain.N+1);
	var.U.row(domain.N+3) = var.U.row(domain.N+2);
}

void MUSCL::initial_conditions(eulerTests &Test){
	//Test.var.X.resize(Test.domain.N+2, 1);
	Test.var.U.resize(Test.domain.N+4, 3);
	Test.var.F.resize(Test.domain.N+2, 3);

	vector eulerL;
	vector eulerR;

	eulerL = Test.var.state_function->conservedVar(Test.initialL);
	eulerR = Test.var.state_function->conservedVar(Test.initialR);

	for (int i=0; i<Test.domain.N; i++){
		Test.domain.X(i+1) = i*Test.domain.dx;
		if (Test.domain.X(i+1)  < Test.x0){
			Test.var.U.row(i+2) = eulerL;
		}
		else Test.var.U.row(i+2) = eulerR;
	}

	boundary_conditions(Test.var, Test.domain);
}

//----------------------------------------------------------------------------
//	Slope Limiters
//----------------------------------------------------------------------------
vector MUSCL::superBee(matrix U, int i){
	//Calculating slope
	vector diMinus = U.row(i) - U.row(i-1);
	vector diPlus = U.row(i+1) - U.row(i);
	vector di = 0.5*diMinus + 0.5*diPlus;

	//calculating limiter
	vector epsiloni(0, 0, 0);
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

/*vector MUSCL::vanLeer(matrix U, int i){

	vector diMinus = U.row(i) - U.row(i-1);
	vector diPlus = U.row(i+1) - U.row(i);
	vector di = 0.5*diMinus + 0.5*diPlus;

	//-----------------------------------------------
	//	Slope limiter -- Van Leer
	//----------------------------------------------

	vector epsiloni(0, 0, 0);
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
}*/

matrix MUSCL::vanLeer(matrix U, int i){

	/*-----------------------------------------------
	 * Slope limiter -- Van Leer
	 ----------------------------------------------*/

	if (U.cols() == 3){

		vector diMinus = U.row(i) - U.row(i-1);
		vector diPlus = U.row(i+1) - U.row(i);
		vector di = 0.5*diMinus + 0.5*diPlus;

		vector epsiloni(0, 0, 0);
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

		matrix diBar(3, 1);
		for (int j=0; j<3; j++){
			 diBar(j) = epsiloni(j)*di(j);
		}

		return diBar;
		//make r = 0 to test first order
	}

	else if (U.cols() == 4){
		vector4 diMinus = U.row(i) - U.row(i-1);
		vector4 diPlus = U.row(i+1) - U.row(i);
		vector4 di = 0.5*diMinus + 0.5*diPlus;

		vector4 epsiloni(0, 0, 0, 0);
		vector4 epsilonR(0, 0, 0, 0);
		vector4 ri(0, 0, 0, 0);

		for (int j=0; j<4; j++){
			ri(j) = diMinus(j)/diPlus(j);
		}

		for (int j=0; j<4; j++){
			if (ri(j) <= 0){
				epsiloni(j) = 0;
			}
			//note: A more refined approach would be to adopt characteristic limiting p510
			else if (ri(j) > 0){
				epsilonR(j) = 2/(1+ri(j));
				epsiloni(j) = fmin(2*ri(j)/(1+ri(j)), epsilonR(j));
			}
		}

		matrix diBar(4, 1);
		for (int j=0; j<4; j++){
			 diBar(j) = epsiloni(j)*di(j);
		}

		return diBar;		
	}

	else {
		throw "Number of conserved variables exceeds 2D";
	}
}

vector MUSCL::minBee(matrix U, int i){


	vector diMinus = U.row(i) - U.row(i-1);
	vector diPlus = U.row(i+1) - U.row(i);
	vector di = 0.5*diMinus + 0.5*diPlus;

	/*-----------------------------------------------
	 * Slope limiter -- MinBee
	 ----------------------------------------------*/

	vector epsiloni(0, 0, 0);
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

void MUSCL::data_reconstruction(matrix U, slopeLimiter a, matrix &ULi, matrix &URi, int N){
	matrix Utmp(0, 0);

	if (ULi.cols() == 3) {
		 Utmp.resize(3, 1);
	}

	if (ULi.cols() == 4) {
		 Utmp.resize(4, 1);
	}
	//slopeLimiter a = getLimiter();

	for (int i=1; i<N+3; i++){ //U goes from 0 to N+3

		//-----------------------------------------------
		// Data Reconstruction
		//----------------------------------------------
		switch(a){
		case MinBee:
			Utmp = U.row(i);
			ULi.row(i) = Utmp - 0.5*minBee(U, i);
			URi.row(i) = Utmp + 0.5*minBee(U, i);
			break;
		case VanLeer:
			Utmp = U.row(i);
			ULi.row(i) = Utmp - (0.5*vanLeer(U, i)).transpose();
			URi.row(i) = Utmp + (0.5*vanLeer(U, i)).transpose();
			break;
		case SuperBee:
			Utmp = U.row(i);
			ULi.row(i) = Utmp - 0.5*superBee(U, i);
			URi.row(i) = Utmp + 0.5*superBee(U, i);
			break;
		case Quit:
			exit(0);
		}

	}
}


void MUSCL::conservative_update_formula(Euler1D &var, Domain1D domain, int i){
	var.U.row(i) = var.U.row(i) - (domain.dt/domain.dx)*(var.F.row(i) - var.F.row(i-1));
}
//F_1 represents F[i-1]
//void MUSCL::conservative_update_formula(vector& U, vector F, vector F_1, double dt, double dx){
//	U = U - (dt/dx)*(F - F_1);
//}

void MUSCL::conservative_update_formula_2D(vector4& U, vector4 F, vector4 F_1, double dt, double dx){
	U = U - (dt/dx)*(F - F_1);
}

//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------

void MUSCL::compute_fluxes(Euler1D &var, Domain1D &domain, int i, matrix ULi, matrix URi, double &Smax){

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
	vector ULbar = ULtmp1 + 0.5*(domain.dt/domain.dx)*(var.state_function->fluxes(ULtmp1) - var.state_function->fluxes(URtmp1)); //UL(i+1)
	vector URbar = URtmp + 0.5*(domain.dt/domain.dx)*(var.state_function->fluxes(ULtmp) - var.state_function->fluxes(URtmp));


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
		Pl = var.state_function->Pressure(hllcUL);
		Pr = var.state_function->Pressure(hllcUR);

		//velocity
		ul = hllcUL(1)/hllcUL(0);
		ur = hllcUR(1)/hllcUR(0);

		//soundspeed
		al = var.state_function->soundspeed(hllcUL);
		ar = var.state_function->soundspeed(hllcUR);


		//Davies Wave Speed Estimates
		double Splus = fmax(abs(ul) + al, abs(ur) + ar);
		SL = -Splus; SR = Splus;
		

				/*
				//---------------------------------------
				// pressure based wave speed estimate
				//---------------------------------------

				//Pressure-based wave speed estimate
				double Ppvrs;
				double Pstar;
				double rhoavg;
				double aavg;
				double ql, qr;

				rhoavg = 0.5*(dl + dr);
				aavg = 0.5*(al + ar);

				Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;

				Pstar = fmax(0.0, Ppvrs);
				if (Pstar <= Pl){
					ql = 1.0;
				}

				else {
					ql = sqrt(1 + ((var.state_function->y+1)/(2*var.state_function->y))*((Pstar/Pl) - 1));
				}

				if (Pstar <= Pr){
					qr = 1.0;
				}

				else {
					qr = sqrt(1 + ((var.state_function->y+1)/(2*var.state_function->y))*((Pstar/Pr) - 1));
				}

				SR = ur + ar*qr;
				SL = ul - al*ql;
				*/


		if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

		Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));
		//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;
		//initialize FL and FR for each timestep
		vector FL(mvl, mvl*ul + Pl, ul*(El + Pl));
		vector FR(mvr, mvr*ur + Pr, ur*(Er + Pr));		

		if (0 <= SL){
			var.F.row(i) = FL;
		}

		//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
		else if (SL<=0 && Sstar>=0){
			double tmpUstar = dl*((SL - ul)/(SL - Sstar));
			vector UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))));
			tmp = var.U.row(i);
			var.F.row(i) = FL + SL*(UstarL - tmp);
		}

		else if (Sstar<=0 && SR>=0){
			double tmpUstar = dr*((SR - ur)/(SR - Sstar));
			vector UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))));
			tmp = var.U.row(i+1);
			var.F.row(i) = FR + SR*(UstarR - tmp);
		}

		else if (0 >= SR){
			var.F.row(i) = FR;
		}
}

void MUSCL::solver(Euler1D &var, Domain1D &domain, double CFL){
	
	matrix ULi(domain.N+4, 3);
	matrix URi(domain.N+4, 3);
	double Smax=0;

	slopeLimiter a = getLimiter();

	double t = 0.0;
	int count = 0;
	do{

		data_reconstruction(var.U, a, ULi, URi, domain.N);
		for (int i=1; i<domain.N+2; i++){
			compute_fluxes(var, domain, i, ULi, URi, Smax);
		}

		//set timestep
		if (count > 5) domain.dt = CFL*(domain.dx/Smax); //updates every timestep
		else domain.dt = 0.2*(domain.dx/Smax);
		if (t + domain.dt > domain.tstop) domain.dt = domain.tstop - t;

		t += domain.dt;
		count += 1;

		//updating U
		for (int i=2; i<domain.N+2; i++){
			conservative_update_formula(var, domain, i);
			//conservative_update_formula(var.U.row(i), var.F.row(i), var.F.row(i-1), domain.dt, domain.dx);
		}

		boundary_conditions(var, domain);

	}while (t < domain.tstop);
	std::cout << count << std::endl;
}

void MUSCL::output(Euler1D &var, Domain1D domain){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	for (int i=2; i<domain.N+2; i++){

		double u = var.U(i, 1)/var.U(i, 0);
		double P = var.state_function->Pressure(var.U, i);
		double e = var.state_function->internalE(var.U, i);

		outfile << domain.X(i-1) << '\t' << var.U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	outfile.close();
	std::cout << "done: MUSCL" << std::endl;
}

void MUSCL::muscl_solver(eulerTests& Test, double CFL){
	initial_conditions(Test);
	solver(Test.var, Test.domain, CFL);
	output(Test.var, Test.domain);
}

//--------------------------------------------------------------------------------
//	2D MUSCL
//--------------------------------------------------------------------------------
//
//
//--------------------------------------------------------------------------------
//	
//--------------------------------------------------------------------------------


void MUSCL::boundary_conditions(Euler2D &var, Domain2D domain){

	if (var.U.cols() == 0 || var.U.rows() == 0){
		throw "Array is empty.";
	}

	//assigning ghost values in the x-direction 
	for (int j=0; j<domain.Ny; j++){
		var.U(1, j+2) = var.U(2, j+2);
		var.U(0, j+2) = var.U(1, j+2);
		var.U(domain.Nx+2, j+2) = var.U(domain.Nx+1, j+2);
		var.U(domain.Nx+3, j+2) = var.U(domain.Nx+2, j+2);
	} 
	//assigning ghost values in the y-direction
	for (int i=0; i<domain.Nx; i++){
		var.U(i+2, 1) = var.U(i+2, 2);
		var.U(i+2, 0) = var.U(i+2, 1);
		var.U(i+2, domain.Ny+2) = var.U(i+2, domain.Ny+1);
		var.U(i+2, domain.Ny+3) = var.U(i+2, domain.Nx+2);
	} 
}

void MUSCL::boundary_conditions_reflective(Euler2D &var, Domain2D domain){

	if (var.U.cols() == 0 || var.U.rows() == 0){
		throw "Array is empty.";
	}

	//assigning ghost values in the x-direction 
	for (int j=0; j<domain.Ny; j++){
		var.U(1, j+2) = var.U(2, j+2);
		var.U(0, j+2) = var.U(1, j+2);
		var.U(domain.Nx+2, j+2) = var.U(domain.Nx+1, j+2);
		var.U(domain.Nx+3, j+2) = var.U(domain.Nx+2, j+2);
		//reflecting the normal velocity
		var.U(1, j+2)(1) = -var.U(2, j+2)(1);
		var.U(0, j+2)(1) = -var.U(1, j+2)(1);
		var.U(domain.Nx+2, j+2)(1) = -var.U(domain.Nx+1, j+2)(1);
		var.U(domain.Nx+3, j+2)(1) = -var.U(domain.Nx+2, j+2)(1);
	} 
	//assigning ghost values in the y-direction
	for (int i=0; i<domain.Nx; i++){
		var.U(i+2, 1) = var.U(i+2, 2);
		var.U(i+2, 0) = var.U(i+2, 1);
		var.U(i+2, domain.Ny+2) = var.U(i+2, domain.Ny+1);
		var.U(i+2, domain.Ny+3) = var.U(i+2, domain.Nx+2);
		//reflecting the normal velocity
		var.U(i+2, 1)(3) = -var.U(i+2, 2)(3);
		var.U(i+2, 0)(3) = -var.U(i+2, 1)(3);
		var.U(i+2, domain.Ny+2)(3) = -var.U(i+2, domain.Ny+1)(3);
		var.U(i+2, domain.Ny+3)(3) = -var.U(i+2, domain.Nx+2)(3);
	} 
}

void MUSCL::initial_conditions(eulerTests2D &Test){
	//Test.var.X.resize(Test.domain.N+2, Test.domain.Ny+2);
	Test.var.U.resize(Test.domain.Nx+4, Test.domain.Ny+4);
	Test.var.F.resize(Test.domain.Nx+2, Test.domain.Ny+2);
	Test.var.G.resize(Test.domain.Nx+2, Test.domain.Ny+2);

	vector4 eulerL;
	vector4 eulerR;

	eulerL = Test.var.state_function->conservedVar2Dx(Test.initialL);
	eulerR = Test.var.state_function->conservedVar2Dx(Test.initialR);

	for (int i=0; i<Test.domain.Nx; i++){
		for (int j=0; j<Test.domain.Ny; j++){
			if (Test.interface(i, j) == false){
				Test.var.U(i+2, j+2) = eulerL;			
			}
			else {
				Test.var.U(i+2, j+2) = eulerR;				
			}
		}
	}
	boundary_conditions(Test.var, Test.domain);
	//Test.var.display(Test.var.U);
}

//Piecewise constant intterpolation -- HLLC, for boundaries with only 1 ghost point
void MUSCL::compute_fluxes(const Euler2D &var, const matrix &U, vector4 &F, int i){
	//Here dx is a proxy for a unit cell in any direction

	//soundspeed
	double al, ar;
	//velocity components
	double ul, ur; //direction currently being swept
	double vl, vr; //other direction

	double dl, dr;
	double Pl, Pr;

	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;

	vector4 tmp;

	/*-------------------------------------------------------
	 * HLLC solver
	 -------------------------------------------------------*/
	vector4 hllcUL = U.row(i);
	vector4 hllcUR = U.row(i+1);

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
		Pl = var.state_function->Pressure(hllcUL);
		Pr = var.state_function->Pressure(hllcUR);

		//velocity
		ul = hllcUL(1)/hllcUL(0);
		ur = hllcUR(1)/hllcUR(0);

		vl = hllcUL(3)/hllcUL(0);
		vr = hllcUR(3)/hllcUR(0);

		//soundspeed
		al = var.state_function->soundspeed(hllcUL);
		ar = var.state_function->soundspeed(hllcUR);


		//Davies Wave Speed Estimates
		//double Splus = fmax(abs(ul) + al, abs(ur) + ar);
		//SL = -Splus; SR = Splus;
		
		double uspecial = (ul + ur)/2. + (al - ar)/(var.state_function->y - 1);
		double aspecial = (al + ar)/2. + (ul - ur)*(var.state_function->y - 1)/4.;
		SL = fmin(0, fmin(ul - al, uspecial - aspecial));
		SR = fmax(0, fmax(ur + ar, uspecial + aspecial));

		Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));
		//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;
		//initialize FL and FR for each timestep
		vector4 FL(mvl, mvl*ul + Pl, ul*(El + Pl), mvl*vl);
		vector4 FR(mvr, mvr*ur + Pr, ur*(Er + Pr), mvr*vr);		

		if (0 <= SL){
			F = FL;
		}

		//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
		else if (SL<=0 && Sstar>=0){
			double tmpUstar = dl*((SL - ul)/(SL - Sstar));
			vector4 UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))),
				tmpUstar*vl);
			tmp = U.row(i);
			F = FL + SL*(UstarL - tmp);
		}

		else if (Sstar<=0 && SR>=0){
			double tmpUstar = dr*((SR - ur)/(SR - Sstar));
			vector4 UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))),
				tmpUstar*vr);
			tmp = U.row(i+1);
			F = FR + SR*(UstarR - tmp);
		}

		else if (0 >= SR){
			F = FR;
		}
}

void MUSCL::compute_fluxes(const Euler2D &var, const Domain2D &domain, const matrix &U, vector4 &F, int i, const matrix &ULi, const matrix &URi, double dx){
	//Here dx is a proxy for a unit cell in any direction

	//soundspeed
	double al, ar;
	//velocity components
	double ul, ur; //direction currently being swept
	double vl, vr; //other direction

	double dl, dr;
	double Pl, Pr;

	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;

	vector4 tmp;

	/*-----------------------------------------------
	 * Evolution by 1/2 time-step
	 ----------------------------------------------*/

	vector4 ULtmp = ULi.row(i);
	vector4 URtmp = URi.row(i);
	vector4 ULtmp1 = ULi.row(i+1);
	vector4 URtmp1 = URi.row(i+1);

	//Post 1/2 time-step evolution left and right states
	vector4 ULbar = ULtmp1 + 0.5*(domain.dt/dx)*(var.state_function->fluxes(ULtmp1) - var.state_function->fluxes(URtmp1)); //UL(i+1)
	vector4 URbar = URtmp + 0.5*(domain.dt/dx)*(var.state_function->fluxes(ULtmp) - var.state_function->fluxes(URtmp));


	/*-------------------------------------------------------
	 * Solution of the piecewise constant Riemann problem pg 180
	 -------------------------------------------------------*/
	/*-------------------------------------------------------
	 * HLLC solver
	 -------------------------------------------------------*/
	vector4 hllcUL = URbar;
	vector4 hllcUR = ULbar;

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
		Pl = var.state_function->Pressure(hllcUL);
		Pr = var.state_function->Pressure(hllcUR);

		//velocity
		ul = hllcUL(1)/hllcUL(0);
		ur = hllcUR(1)/hllcUR(0);

		vl = hllcUL(3)/hllcUL(0);
		vr = hllcUR(3)/hllcUR(0);

		//soundspeed
		al = var.state_function->soundspeed(hllcUL);
		ar = var.state_function->soundspeed(hllcUR);


		//Davies Wave Speed Estimates
		//double Splus = fmax(abs(ul) + al, abs(ur) + ar);
		//SL = -Splus; SR = Splus;
		
		double uspecial = (ul + ur)/2. + (al - ar)/(var.state_function->y - 1);
		double aspecial = (al + ar)/2. + (ul - ur)*(var.state_function->y - 1)/4.;
		SL = fmin(0, fmin(ul - al, uspecial - aspecial));
		SR = fmax(0, fmax(ur + ar, uspecial + aspecial));

				
				//---------------------------------------
				// pressure based wave speed estimate
				//---------------------------------------
				/*
				//Pressure-based wave speed estimate
				double Ppvrs;
				double Pstar;
				double rhoavg;
				double aavg;
				double ql, qr;

				rhoavg = 0.5*(dl + dr);
				aavg = 0.5*(al + ar);

				Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;

				Pstar = fmax(0.0, Ppvrs);
				if (Pstar <= Pl){
					ql = 1.0;
				}

				else {
					ql = sqrt(1 + ((var.state_function->y+1)/(2*var.state_function->y))*((Pstar/Pl) - 1));
				}

				if (Pstar <= Pr){
					qr = 1.0;
				}

				else {
					qr = sqrt(1 + ((var.state_function->y+1)/(2*var.state_function->y))*((Pstar/Pr) - 1));
				}

				SR = ur + ar*qr;
				SL = ul - al*ql;
				*/

		//if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

		Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));
		//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;
		//initialize FL and FR for each timestep
		vector4 FL(mvl, mvl*ul + Pl, ul*(El + Pl), mvl*vl);
		vector4 FR(mvr, mvr*ur + Pr, ur*(Er + Pr), mvr*vr);		

		if (0 <= SL){
			F = FL;
		}

		//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
		else if (SL<=0 && Sstar>=0){
			double tmpUstar = dl*((SL - ul)/(SL - Sstar));
			vector4 UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))),
				tmpUstar*vl);
			tmp = U.row(i);
			F = FL + SL*(UstarL - tmp);
		}

		else if (Sstar<=0 && SR>=0){
			double tmpUstar = dr*((SR - ur)/(SR - Sstar));
			vector4 UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))),
				tmpUstar*vr);
			tmp = U.row(i+1);
			F = FR + SR*(UstarR - tmp);
		}

		else if (0 >= SR){
			F = FR;
		}
}


void MUSCL::solver(Euler2D &var, Domain2D &domain, double CFL){
	
	matrix ULix(domain.Nx+4, 4);
	matrix URix(domain.Nx+4, 4);

	matrix ULiy(domain.Ny+4, 4);
	matrix URiy(domain.Ny+4, 4);

	//Temporary storage for x and y conserved variables 
	matrix Uxn(domain.Nx+4, 4);
	matrix Uyn(domain.Ny+4, 4);


	slopeLimiter a = VanLeer;//getLimiter();

	double t = 0.0;
	int count = 0;
	do{
		//set timestep, taking maximum wavespeed in x or y directions
		double Smax=0;
		//double Smax_y=0;
		for (int i=2; i<domain.Nx+2; i++){
			for (int j=2; j<domain.Ny+2; j++){
				vector4 Utmp = var.U(i, j);
				double a_ij = var.state_function->soundspeed(Utmp);
				double d_ij_x = Utmp(1)/Utmp(0); //velocity in x direction
				double d_ij_y = Utmp(3)/Utmp(0); //velocity in y direction
				double vn = sqrt(pow(d_ij_x, 2) + pow(d_ij_y, 2));
				double lambda = fmax(fmax(abs(d_ij_x) + a_ij, abs(d_ij_y) + a_ij), vn + a_ij);
				if (lambda > Smax) Smax = lambda;
			}
		}
		if (count <= 5){
				//using a CFL number of 0,.2 for the first 5 timesteps
				domain.dt = 0.2*fmin((domain.dx/Smax), (domain.dy/Smax));
			}
		else {
				domain.dt = CFL*fmin((domain.dx/Smax), (domain.dy/Smax));
		}
		//domain.dt = CFL*fmin((domain.dx/Smax_xtest), (domain.dy/Smax_ytest));
		if (t + domain.dt > domain.tstop) domain.dt = domain.tstop - t;

		t += domain.dt;
		count += 1;
		//sweeping in the x-direction for each y row
		for (int j=0; j<domain.Ny; j++){
			//Within each row, store the values for each grid point in x with a single row list
			for(int i=0; i<domain.Nx+4; i++){
				Uxn.row(i) = var.U(i, j+2);
				//the sweep in x is only performed within the real y grid
			} 
			//calculate and update to new interpolated values U_i
			data_reconstruction(Uxn, a, ULix, URix, domain.Nx);

			//sweeping in the x-direction for each y row
			for (int i=1; i<domain.Nx+2; i++){
				vector4 Fx(0, 0, 0, 0);
				compute_fluxes(var, domain, Uxn, Fx, i, ULix, URix, domain.dx); //compute flux needs to change
				var.F(i, j+1) = Fx; //storing the computed flux.
			}
		}
		//sweeping in the y-direction for each x row
		for (int i=0; i<domain.Nx; i++){
			for (int j=0; j<domain.Ny+4; j++){
				Uyn.row(j) = var.swap_xy(var.U(i+2, j));
			}
			data_reconstruction(Uyn, a, ULiy, URiy, domain.Ny);

			for (int j=1; j<domain.Ny+2; j++){
				vector4 Fy(0, 0, 0, 0);
				compute_fluxes(var, domain, Uyn, Fy, j, ULiy, URiy, domain.dy); //compute flux needs to change
				var.G(i+1, j) = Fy; //storing the computed flux.		
			}
		}

		//updating U
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				conservative_update_formula_2D(var.U(i, j+2), var.F(i, j+1), var.F(i-1, j+1), domain.dt/2, domain.dx);
			}
		}
		for (int i=0; i<domain.Nx; i++){
			for (int j=2; j<domain.Ny+2; j++){
				conservative_update_formula_2D(var.U(i+2, j), var.swap_xy(var.G(i+1, j)), var.swap_xy(var.G(i+1, j-1)), domain.dt, domain.dy);
			}
		}		
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				conservative_update_formula_2D(var.U(i, j+2), var.F(i, j+1), var.F(i-1, j+1), domain.dt/2, domain.dx);
			}
		}
		boundary_conditions(var, domain);

		if (count%50==0) std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << Smax_x << '\t' << Smax_y << std::endl;
		//if (count == 2) t = domain.tstop;
	}while (t < domain.tstop);
	std::cout << count << std::endl;
}

void MUSCL::output(const Euler2D &var, const Domain2D &domain){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	for (int i=2; i<domain.Nx+2; i++){
		for (int j=2; j<domain.Ny+2; j++){
			vector4 Ux = var.U(i, j);
			double u = sqrt(pow(Ux(1)/Ux(0), 2) + pow(Ux(3)/Ux(0), 2));
			double P = var.state_function->Pressure(Ux);
			double e = var.state_function->internalE(Ux);

			//central difference to calculate partial derivatives in x, y for density
			double grad_density_x = (var.U(i+1, j)(0) - var.U(i-1, j)(0))/(2*domain.dx);
			double grad_density_y = (var.U(i, j+1)(0) - var.U(i, j-1)(0))/(2*domain.dy);
			//calculating the numerical schlieren
			double schlieren = exp((-20*sqrt(pow(grad_density_x, 2) + pow(grad_density_y, 2)))/(1000*Ux(0)));

			outfile << domain.dx*(i-2) << '\t' << domain.dy*(j-2) << '\t' << Ux(0) << '\t' << u
			<< '\t' << P << '\t' << e << '\t' << schlieren << std::endl;
		}
		outfile << std::endl;
	}
	//plotting a slice
	/*for (int j=(domain.Ny/2 +2); j<domain.Ny+2; j++){
		vector4 Ux = var.U(domain.Nx/2, j);
		double u = sqrt(pow(Ux(1)/Ux(0), 2) + pow(Ux(3)/Ux(0), 2));
		double P = var.state_function->Pressure(Ux);
		double e = var.state_function->internalE(Ux);

		outfile << domain.dy*(j-2) << '\t' << Ux(0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}*/
	outfile.close();
	std::cout << "done: MUSCL" << std::endl;
}

void MUSCL::muscl_solver(eulerTests2D& Test, double CFL){
	initial_conditions(Test);
	solver(Test.var, Test.domain, CFL);
	output(Test.var, Test.domain);
}


//--------------------------------------------------------------------------------
//	Exact Solver
//--------------------------------------------------------------------------------
//Note, W contains primitive variables (density, velocity, pressure)

/*EXACT::EXACT(int N, double dx, double x0, double y, vector WL, vector WR)
	: dx(dx), x0(x0), W(N, 3), WL(WL), WR(WR), TOL(1e-6), y(y), cL(0), cR(0), 
	CONST1(0),CONST2(0), CONST3(0), CONST4(0), CONST5(0), CONST6(0), CONST7(0), CONST8(0) {}*/

void EXACT::initial_conditions(eulerTests& Test){
	WL = Test.initialL;
	WR = Test.initialR;

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

void EXACT::solver(eulerTests &Test){
	initial_conditions(Test);
	sampling(Test.domain.tstop);
	output();
}









