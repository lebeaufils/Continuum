#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

const int N = 100; // number of cells
const double L = 1.0;
const double x0 = 0.5;
const double dx = L/N; //length of domain
const double CFL = 0.9;
const double y = 1.4; //heat capacity ratio Cp/Cv

const Eigen::Vector4d initialL(1.0, 2.0, 1.0, 1.0); //density, velocity, pressure
const Eigen::Vector4d initialR(1.0, 2.0, 1.0, 0.0); //left and right states

typedef Eigen::Matrix<double, N+2, 4> MatrixN2;
typedef Eigen::Matrix<double, N+1, 4> MatrixN1;
typedef Eigen::Matrix<double, N+2, 1> VectorN2;
typedef Eigen::Matrix<double, N+1, 1> VectorN1;

void boundary_conditions(MatrixN2 &U){
	U.row(0) = U.row(1);
	U.row(N+1) = U.row(N);
}

void initial_conditions(VectorN2 &X, MatrixN2 &U){

	Eigen::Vector4d eulerL;
	Eigen::Vector4d eulerR;

	//rho
	eulerL(0) = initialL(0);
	eulerR(0) = initialR(0);

	//rhou
	eulerL(1) = initialL(0)*initialL(1);
	eulerR(1) = initialR(0)*initialR(1);

	//E
	eulerL(2) = initialL(2)/(y-1) + 0.5*initialL(0)*initialL(1)*initialL(1);
	eulerR(2) = initialR(2)/(y-1) + 0.5*initialR(0)*initialR(1)*initialR(1);

	//advected tracer
	eulerL(3) = initialL(3);
	eulerR(3) = initialR(3);

	for (int i=0; i<N+1; i++){
		X(i) = i*dx;
		if (X(i)  < x0){
			U.row(i+1) = eulerL;
		}
		else U.row(i+1) = eulerR;
	}

	boundary_conditions(U);
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
	//augmented addition
	double lambdal, lambdar;


	double mvl, mvr;
	double El, Er;

	double SL, SR, Sstar;

	double Smax = 0;
	MatrixN1 F;
	Eigen::Vector4d tmp;

	//MatrixN2 U_1;

	double dt = 0.01;
	double t = 0.0;
	int count = 0;
	do{

		for (int i=0; i<N+1; i++){

			//if (count == 1) std::cout << U.row(i)<< std::endl;
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
			Pl = (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0));
			Pr = (y-1)*(U(i+1, 2) - 0.5*U(i+1, 0)*pow((U(i+1, 1)/U(i+1, 0)),2.0));

			//velocity
			ul = U(i, 1)/U(i, 0);
			ur = U(i+1, 1)/U(i+1, 0);

			//lambda
			lambdal = U(i, 3)/U(i, 0);
			lambdar = U(i+1, 3)/U(i+1, 0);

			//soundspeed
			al = sqrt(y*(Pl/dl));
			ar = sqrt(y*(Pr/dr));

			/*
			 * pressure based wave speed estimate
			 */

			rhoavg = 0.5*(dl + dr);
			aavg = 0.5*(al + ar);

			Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;

			std::cout << lambdal << '\t' << lambdar << '\t' << U(i, 3) << std::endl;

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
			//if (count == 1) std::cout << Pr << '\t' << Pl  << '\t' << dr << '\t' << dl<< std::endl;

			//initialize FL and FR for each timestep
			Eigen::Vector4d FL(mvl, mvl*ul + Pl, ul*(El + Pl), mvl*lambdal);
			Eigen::Vector4d FR(mvr, mvr*ur + Pr, ur*(Er + Pr), mvr*lambdar);

			if (0 <= SL){
				F.row(i) = FL;
			}

			//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
			else if (SL<=0 && Sstar>=0){
				double tmpUstar = dl*((SL - ul)/(SL - Sstar));
				Eigen::Vector4d UstarL(tmpUstar, tmpUstar*Sstar, tmpUstar*((El/dl) + (Sstar - ul)*(Sstar + (Pl/(dl*(SL - ul))))), tmpUstar*lambdal);
				tmp = U.row(i);
				F.row(i) = FL + SL*(UstarL - tmp);
			}

			else if (Sstar<=0 && SR>=0){
				double tmpUstar = dr*((SR - ur)/(SR - Sstar));
				Eigen::Vector4d UstarR(tmpUstar, tmpUstar*Sstar, tmpUstar*((Er/dr) + (Sstar - ur)*(Sstar + (Pr/(dr*(SR - ur))))), tmpUstar*lambdar);
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
		double P = (y-1)*(U(i+1, 2) - 0.5*U(i+1, 0)*pow((U(i+1, 1)/U(i+1, 0)),2.0));
		//double e = (y-1)*(U(i+1, 2)-0.5*U(i+1, 1)*(U(i+1, 1)/U(i+1, 0)))/(U(i,0)*(y-1));
		//double e = P/(U(i+1, 0)*(y - 1));
		double lambda = U(i+1, 3)/U(i+1, 0);

		outfile << X(i) << '\t' << U(i+1, 0) << '\t' << u
				<< '\t' << P << '\t' << lambda << std::endl;
	}
	std::cout << "done" << std::endl;
}

int main(void){

	VectorN2 X; //mesh
	MatrixN2 U; //U(i)(0) = rho | U(i)(1) = rhou | U(i)(2) = E

	initial_conditions(X, U);
	double tstop = 0.1;
	HLLCsolver(U, tstop, CFL);
	output(U, X, 0.2);
}




