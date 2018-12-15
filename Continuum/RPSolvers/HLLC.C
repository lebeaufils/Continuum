/*
 * HLLC.C
 *
 *  Created on: 9 Dec 2018
 *      Author: forte
 */

#include "HLLC.h"


HLLC::HLLC(double c, eulerTests Test)
	: CFL(c), X(Test.N+2, 1), U(Test.N+2, 3), F(Test.N+1, 3), count(0), N(Test.N), dx(Test.L/Test.N),dt(0){
}

void HLLC::boundary_conditions(){
	U.row(0) = U.row(1);
	U.row(N+1) = U.row(N);
}

void HLLC::initial_conditions(eulerTests Test, IdealGas IG){
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

void HLLC::solver(IdealGas IG, eulerTests Test){

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
	double Smax = 0;

	vector tmp;

	double t = 0.0;
	do{
		for (int i=0; i<N+1; i++){

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
			Pl = IG.Pressure(U, i);
			Pr = IG.Pressure(U, i+1);

			//velocity
			ul = U(i, 1)/U(i, 0);
			ur = U(i+1, 1)/U(i+1, 0);

			//soundspeed
			al = IG.soundspeed(U, i);
			ar = IG.soundspeed(U, i+1);

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
	std::cout << t << std::endl;
	std::cout << count << std::endl;
}

void HLLC::output(IdealGas IG){

	std::ofstream outfile;
	outfile.open("dataHLLC.txt");

	for (int i=1; i<N+2; i++){

		double u = U(i, 1)/U(i, 0);
		double P = IG.Pressure(U, i);
		double e = IG.internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}

/*----------------------------------------------------------------------------
 * JWL
 ----------------------------------------------------------------------------*/

void HLLC::initial_conditions(eulerTests Test, JWL MG){
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

void HLLC::solver(JWL MG, eulerTests Test){

	double al, ar;

	double ul, ur; //note u is the eigenvalue of the euler equation
	double dl, dr;
	double Pl, Pr;
	
	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;
	double Splus, Smax = 0;

	vector tmp;

	double t = 0.0;
	do{

		for (int i=0; i<N+1; i++){
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

			//derived variables
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
			SL = -Splus;
			SR = Splus;

			//finding Smax for the whole domain (for each timestep)
			if (Splus > Smax) Smax = Splus;

			Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));

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

	}while(t < Test.tstop);
	std::cout << count << std::endl;
}

void HLLC::output(JWL MG){

	std::ofstream outfile;
	outfile.open("dataHLLC.txt");

	for (int i=1; i<N+2; i++){

		double u = U(i, 1)/U(i, 0);
		double P = MG.Pressure(U, i);
		double e = MG.internalE(U, i);

		outfile << X(i) << '\t' << U(i, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}


