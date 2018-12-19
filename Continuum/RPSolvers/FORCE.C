/*
 * FORCE.C
 *
 *  Created on: 19 Dec 2018
 *      Author: forte
 */
#include "FORCE.h"

FORCE::FORCE(double c, eulerTests Test)
	:CFL(c), X(Test.N+2, 1), U(Test.N+2, 3), F(Test.N+1, 3), count(0), N(Test.N), dx(Test.L/Test.N), dt(0){
}

void FORCE::boundary_conditions(){
	U.row(0) = U.row(1);
	U.row(N+1) = U.row(N);
}

void FORCE::initial_conditions(IdealGas IG, eulerTests Test){
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

	for (int i=0; i<N+1; i++){
		X(i+1) = i*dx;
		if (X(i+1)  < Test.x0){
			U.row(i+1) = eulerL;
		}
		else U.row(i+1) = eulerR;
	}

	boundary_conditions();
}

void FORCE::initial_conditions(JWL MG, eulerTests Test){
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

	for (int i=0; i<N+1; i++){
		X(i+1) = i*dx;
		if (X(i+1)  < Test.x0){
			U.row(i+1) = eulerL;
		}
		else U.row(i+1) = eulerR;
	}
	boundary_conditions();
}

void FORCE::solver(IdealGas IG, eulerTests Test){

	double al, ar; //soundspeed
	double ul, ur; //note u is the eigenvalue of the euler equation
	double dl, dr;
	double Pl, Pr;
	double mvl, mvr;
	double El, Er;

	double Splus, Smax = 0; //max soundspeed

	vector tmpl, tmpr;

	double t = 0.0;
	do{

		for (int i=0; i<(N+1); i++){

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

			//Pressure
			Pl = IG.Pressure(U, i);
			Pr = IG.Pressure(U, i+1);

			//velocity
			ul = U(i, 1)/U(i, 0);
			ur = U(i+1, 1)/U(i+1, 0);
			//soundspeed
			al = IG.soundspeed(U, i);
			ar = IG.soundspeed(U, i+1);

			/*--------------------------------------------------------------
			 * Davies wave speed estimate
			 *--------------------------------------------------------------*/

			//S+ = max{| uL | +aL,| uR | +aR} .
			Splus = std::max(abs(ul) + al, abs(ur) + ar);

			//finding Smax for the whole domain (for each timestep)
			if (Splus > Smax) Smax = Splus;

			/*-----------------------------------------------
			 * Important for LxF Flux
			 -----------------------------------------------*/
			//ensure dt is non zero for first time step
			if (count == 0) dt = CFL*(dx/Smax);

			//updating flux
			vector FUl(mvl, mvl*ul+Pl, ul*(El+Pl));
			vector FUr(mvr, mvr*ur+Pr, ur*(Er+Pr));

			tmpl = U.row(i); tmpr = U.row(i+1);

			vector LxF = 0.5*(FUl + FUr) + 0.5*(dx/dt)*(tmpl - tmpr);
			vector R = 0.5*(tmpl + tmpr) + 0.5*(dt/dx)*(FUl - FUr);
			vector FR = IG.f(R);

			F.row(i) = 0.5*(LxF + FR);
		}
		//end of domain loop
		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > Test.tstop) dt = Test.tstop - t;

		//updating U
		for (int i=1; i<N+1; i++){
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
		}
		boundary_conditions();

		t += dt;
		count += 1;
	}while (t < Test.tstop);
	std::cout << count << std::endl;
}

void FORCE::solver(JWL MG, eulerTests Test){

	double al, ar; //soundspeed
	double ul, ur; //note u is the eigenvalue of the euler equation
	double dl, dr;
	double Pl, Pr;
	double mvl, mvr;
	double El, Er;

	double Splus, Smax = 0; //max soundspeed

	vector tmpl, tmpr;

	double t = 0.0;
	do{

		for (int i=0; i<(N+1); i++){

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

			//Pressure
			Pl = MG.Pressure(U, i);
			Pr = MG.Pressure(U, i+1);

			//velocity
			ul = U(i, 1)/U(i, 0);
			ur = U(i+1, 1)/U(i+1, 0);
			//soundspeed
			al = MG.soundspeed(U, i);
			ar = MG.soundspeed(U, i+1);

			/*--------------------------------------------------------------
			 * Davies wave speed estimate
			 *--------------------------------------------------------------*/

			//S+ = max{| uL | +aL,| uR | +aR} .
			Splus = std::max(abs(ul) + al, abs(ur) + ar);

			//finding Smax for the whole domain (for each timestep)
			if (Splus > Smax) Smax = Splus;

			/*-----------------------------------------------
			 * Important for LxF Flux
			 -----------------------------------------------*/
			//ensure dt is non zero for first time step
			if (count == 0) dt = CFL*(dx/Smax);

			//updating flux
			vector FUl(mvl, mvl*ul+Pl, ul*(El+Pl));
			vector FUr(mvr, mvr*ur+Pr, ur*(Er+Pr));

			tmpl = U.row(i); tmpr = U.row(i+1);

			vector LxF = 0.5*(FUl + FUr) + 0.5*(dx/dt)*(tmpl - tmpr);
			vector R = 0.5*(tmpl + tmpr) + 0.5*(dt/dx)*(FUl - FUr);
			vector FR = MG.f(R);

			F.row(i) = 0.5*(LxF + FR);
		}
		//end of domain loop
		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > Test.tstop) dt = Test.tstop - t;

		//updating U
		for (int i=1; i<N+1; i++){
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
		}
		boundary_conditions();

		t += dt;
		count += 1;
	}while (t < Test.tstop);
	std::cout << count << std::endl;
}

void FORCE::output(IdealGas IG){

	std::ofstream outfile;
	outfile.open("dataFORCE.txt");

	for (int i=1; i<N+2; i++){

		double u = U(i, 1)/U(i, 0);
		double P = IG.Pressure(U, i);
		double e = IG.internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}

void FORCE::output(JWL MG){

	std::ofstream outfile;
	outfile.open("dataFORCE.txt");

	for (int i=1; i<N+2; i++){

		double u = U(i, 1)/U(i, 0);
		double P = MG.Pressure(U, i);
		double e = MG.internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}


