#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

const int N = 100; // number of cells
const double L = 1.0; //1m
const double x0 = 0.5;
const double dx = L/N; //length of domain
const double CFL = 0.9;
const double y = 1.4; //heat capacity ratio Cp/Cv
const double tstop = 12e-6;

//Mie-Gruneisen constants
const double Gruneisen = 0.25;
const double d0 = 1840; //kg m-3
const double A =  854.5e9; //GPa
const double B = 20.5e9; //GPa
const double R1 = 4.6;
const double R2 = 1.35;

const Eigen::Vector3d initialL(1700, 0.0, 1e12); //density, velocity, pressure
const Eigen::Vector3d initialR(1000, 0.0, 5e10); //left and right states
//const Eigen::Vector3d initialL(1.0, 0.75, 1.0); //density, velocity, pressure
//const Eigen::Vector3d initialR(0.125, 0.0, 0.1); //left and right states


typedef Eigen::Matrix<double, N+2, 3> MatrixN2;
typedef Eigen::Matrix<double, N+1, 3> MatrixN1;
typedef Eigen::Matrix<double, N+2, 1> VectorN2;
typedef Eigen::Matrix<double, N+1, 1> VectorN1;

void boundary_conditions(MatrixN2 &U){
	U.row(0) = U.row(1);
	U.row(N+1) = U.row(N);
}

void initial_conditions(VectorN2 &X, MatrixN2 &U){

	Eigen::Vector3d eulerL;
	Eigen::Vector3d eulerR;

	//rho
	eulerL(0) = initialL(0);
	eulerR(0) = initialR(0);

	//rhou
	eulerL(1) = initialL(0)*initialL(1);
	eulerR(1) = initialR(0)*initialR(1);

	//E
	eulerL(2) = initialL(0)*(0.5*pow(initialL(1), 2.0) + initialL(2)/(initialL(0)*Gruneisen));
	eulerR(2) = initialR(0)*(0.5*pow(initialR(1), 2.0) + initialR(2)/(initialR(0)*Gruneisen));

	for (int i=0; i<N+2; i++){
		X(i) = i*dx;
		if (X(i)  < x0){
			U.row(i) = eulerL;
		}
		else U.row(i) = eulerR;
	}

	boundary_conditions(U);
}

double Pref(double d){ //function of density
	double dprime = d0/d;

	double Pref = A*exp(-R1*dprime) + B*exp(-R2*dprime);
	return Pref;
}

double eref(double d){
	double dprime = d0/d;

	double eref = A/(R1*d0)*exp(-R1*dprime) + B/(R2*d0)*exp(-R2*dprime);
	return eref;
}

double internalE(MatrixN2 U, int i){
	//E = ρ ( 0.5*V2 + e)
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double Pressure(MatrixN2 U, int i){
	double Pressure = Pref(U(i, 0)) + Gruneisen*U(i, 0)*(internalE(U, i) - eref(U(i, 0)));
	return Pressure;
}

double soundspeed(MatrixN2 U, int i){
	double d = U(i, 0);
	double dprime = d0/d;
	double e = internalE(U, i);
	double csquare = (e - eref(d))*Gruneisen*(Gruneisen + 1)
			+ (1/pow(d, 2.0))*(A*R1*d0*pow(e, -R1*dprime) - B*R2*d0*pow(e, -R2*dprime));
	return sqrt(csquare);
}

void HLLCsolver(MatrixN2 &U, double tstop, double CFL){

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
	MatrixN1 F;
	Eigen::Vector3d tmp;

	//MatrixN2 U_1;

	//JWLparameters

	double dt = 0;
	double t = 0.0;
	int count = 0;
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
			Pl = Pressure(U, i);
			Pr = Pressure(U, i+1);

			//velocity
			ul = U(i, 1)/U(i, 0);
			ur = U(i+1, 1)/U(i+1, 0);

			//soundspeed
			al = soundspeed(U, i);
			ar = soundspeed(U, i+1);

			/*
			 * pressure based wave speed estimate
			 */

			rhoavg = 0.5*(dl + dr);
			aavg = 0.5*(al + ar);

			Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;

			//std::cout << rhoavg << '\t' << aavg << '\t' << Ppvrs << '\t' << dr << std::endl;

			//if (Ppvrs >= 0) Pstar = Ppvrs;
			//else if (Ppvrs < 0) Pstar = 0;
			Pstar = fmax(0.0, Ppvrs);

			if (Pstar <= Pl){
				ql = 1.0;
			}

			else if (Pstar > Pl){
				ql = sqrt(1 + ((y+1)/(2*y))*((Pstar/Pl) - 1));
			}

			if (Pstar <= Pr){
				qr = 1.0;
			}

			else if (Pstar > Pr){
				qr = sqrt(1 + ((y+1)/(2*y))*((Pstar/Pr) - 1));
			}

			SR = ur + ar*qr;
			SL = ul - al*ql;
			//finding Smax for the whole domain (for each timestep)
			if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

			Sstar = (Pr - Pl + dl*ul*(SL - ul) - dr*ur*(SR - ur))/(dl*(SL - ul) - dr*(SR - ur));

			//initialize FL and FR for each timestep
			Eigen::Vector3d FL(mvl, mvl*ul + Pl, ul*(El + Pl));
			Eigen::Vector3d FR(mvr, mvr*ur + Pr, ur*(Er + Pr));

			if (0 <= SL){
				F.row(i) = FL;
			}

			//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
			else if (SL<=0 && Sstar>=0){
				double tmpUstar = dl*((SL - ul)/(SL - Sstar));
				Eigen::Vector3d UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))));
				tmp = U.row(i);
				F.row(i) = FL + SL*(UstarL - tmp);
			}

			else if (Sstar<=0 && SR>=0){
				double tmpUstar = dr*((SR - ur)/(SR - Sstar));
				Eigen::Vector3d UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))));
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
		if (t + dt > tstop) dt = tstop - t;

		//updating U
		for (int i=1; i<N+1; i++){
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
		}
		boundary_conditions(U);

		t += dt;
		count += 1;

	}while(t < tstop);
	std::cout << count << std::endl;
}

void output(MatrixN2 U, VectorN2 X, double tstop){

	std::stringstream filename;
	filename << "Euler_" << std::setw(3) << std::setfill('0') << tstop << "s.txt";

	std::ofstream outfile;
	outfile.open(filename.str().c_str());

	for (int i=0; i<(N+1); i++){

		double u = U(i+1, 1)/U(i+1, 0);
		//double P = (y-1)*(U(i+1, 2) - 0.5*U(i+1, 0)*pow((U(i+1, 1)/U(i+1, 0)),2.0));
		double P = Pressure(U, i+1);
		//double e = P/(U(i+1, 0)*(y - 1));
		//double e = internalE(U, i+1);

		double c = soundspeed(U, i+1);

		outfile << X(i) << '\t' << U(i+1, 0) << '\t' << u
				<< '\t' << P << '\t' << c << std::endl;
	}
	std::cout << "done" << std::endl;
}

int main(void){

	VectorN2 X; //mesh
	MatrixN2 U; //U(i)(0) = rho | U(i)(1) = rhou | U(i)(2) = E

	initial_conditions(X, U);
	HLLCsolver(U, tstop, CFL);
	output(U, X, 0.2);
}




