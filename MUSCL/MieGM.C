/*
 * MUSCL hancock method
 */

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
const double tstop = 12e-6;

//Mie-Gruneisen constants
const double Gruneisen = 0.25;
const double d0 = 1840; //kg m-3
const double A =  854.5e9; //GPa
const double B = 20.5e9; //GPa
const double R1 = 4.6;
const double R2 = 1.35;

//initial properties of explosive
const Eigen::Vector3d initialL(1700, 0.0, 1e12); //density, velocity, pressure
const Eigen::Vector3d initialR(1000, 0.0, 5e10); //left and right states


typedef Eigen::Matrix<double, N+4, 3> MatrixN3;
typedef Eigen::Matrix<double, N+2, 3> MatrixN1;
typedef Eigen::Matrix<double, N+2, 1> VectorN2;
//typedef Eigen::Matrix<double, N+1, 1> VectorN1;

void boundary_conditions(MatrixN3 &U){
	U.row(1) = U.row(2); //muscl requires an extra ghost cell
	U.row(0) = U.row(1);
	U.row(N+2) = U.row(N+1);
	U.row(N+3) = U.row(N+2);
}

void initial_conditions(VectorN2 &X, MatrixN3 &U){

	Eigen::Vector3d eulerL;
	Eigen::Vector3d eulerR;

	//rho
	eulerL(0) = initialL(0);
	eulerR(0) = initialR(0);

	//rhou
	eulerL(1) = initialL(0)*initialL(1);
	eulerR(1) = initialR(0)*initialR(1);

	//M-G version of E.. please check
	//E
	eulerL(2) = initialL(0)*(0.5*pow(initialL(1), 2.0) + initialL(2)/(initialL(0)*Gruneisen));
	eulerR(2) = initialR(0)*(0.5*pow(initialR(1), 2.0) + initialR(2)/(initialR(0)*Gruneisen));

	for (int i=0; i<N+2; i++){
		X(i) = i*dx;
		if (X(i)  < x0){
			U.row(i+1) = eulerL;
		}
		else U.row(i+1) = eulerR;
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

double internalE(MatrixN3 U, int i){
	//E = ρ ( 0.5*V2 + e)
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double Pressure(MatrixN3 U, int i){
	double Pressure = Pref(U(i, 0)) + Gruneisen*U(i, 0)*(internalE(U, i) - eref(U(i, 0)));
	return Pressure;
}

double PressureScalar(Eigen::Vector3d U){
	double e = U(2)/U(0) - 0.5*pow(U(1)/U(0), 2.0);
	double Pressure = Pref(U(0)) + Gruneisen*U(0)*(e - eref(U(0)));
	return Pressure;
}

double soundspeed(MatrixN3 U, int i){
	double d = U(i, 0);
	double dprime = d0/d;
	double e = internalE(U, i);
	double csquare = (e - eref(d))*Gruneisen*(Gruneisen + 1)
			+ (1/pow(d, 2.0))*(A*R1*d0*pow(e, -R1*dprime) - B*R2*d0*pow(e, -R2*dprime));
	return sqrt(csquare);
}

double soundspeedScalar(Eigen::Vector3d U){
	double d = U(0);
	double dprime = d0/d;
	double e = U(2)/U(0) - 0.5*pow(U(1)/U(0), 2.0);
	double csquare = (e - eref(d))*Gruneisen*(Gruneisen + 1)
			+ (1/pow(d, 2.0))*(A*R1*d0*pow(e, -R1*dprime) - B*R2*d0*pow(e, -R2*dprime));
	return sqrt(csquare);
}


Eigen::Vector3d f(Eigen::Vector3d U){
	Eigen::Vector3d flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + PressureScalar(U);
	flux(2) = (U(1)/U(0))*(U(2) + PressureScalar(U));
	return flux;
}


void MHM(MatrixN3 &U, double tstop, double CFL){

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

	//MUSCL
	MatrixN3 ULi;
	MatrixN3 URi;

	MatrixN1 Ftmp;
	MatrixN3 Utmp;

	double dt = 0.01;
	double t = 0.0;
	int count = 0;
	do{

		for (int i=1; i<N+3; i++){ //U goes from 0 to N+3

			Eigen::Vector3d diMinus = U.row(i) - U.row(i-1);
			Eigen::Vector3d diPlus = U.row(i+1) - U.row(i);
			Eigen::Vector3d di = 0.5*diMinus + 0.5*diPlus;

			/*-----------------------------------------------
			 * Slope limiter -- SuperBee
			 ----------------------------------------------*/
/*
			Eigen::Vector3d epsiloni(1, 1, 1);
			Eigen::Vector3d epsilonR(0, 0, 0);
			Eigen::Vector3d ri(0, 0, 0);

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

			Eigen::Vector3d diBar(0, 0, 0);
			for (int j=0; j<3; j++){
				 diBar(j) = epsiloni(j)*di(j);
			}
			//std::cout << diBar.transpose() << std::endl;

*/
			/*-----------------------------------------------
			 * Slope limiter -- Van Leer
			 ----------------------------------------------*/

			Eigen::Vector3d epsiloni(0, 0, 0);
			Eigen::Vector3d epsilonR(0, 0, 0);
			Eigen::Vector3d ri(0, 0, 0);

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

			Eigen::Vector3d diBar(0, 0, 0);
			for (int j=0; j<3; j++){
				 diBar(j) = epsiloni(j)*di(j);
			}


			//make r = 0 to test first order
			/*-----------------------------------------------
			 * Slope limiter -- MinBee
			 ----------------------------------------------*/
			/*
			Eigen::Vector3d epsiloni(0, 0, 0);
			Eigen::Vector3d epsilonR(0, 0, 0);
			Eigen::Vector3d ri(0, 0, 0);

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

			Eigen::Vector3d diBar(0, 0, 0);
			for (int j=0; j<3; j++){
				 diBar(j) = epsiloni(j)*di(j);
			}
			*/

			/*-----------------------------------------------
			 * Data Reconstruction
			 ----------------------------------------------*/

			Eigen::Vector3d Utmp = U.row(i);
			ULi.row(i) = Utmp - 0.5*diBar;
			URi.row(i) = Utmp + 0.5*diBar;
		}

			/*-----------------------------------------------
			 * Evolution by 1/2 time-step
			 ----------------------------------------------*/
		for (int i=1; i<N+2; i++){

			Eigen::Vector3d ULtmp = ULi.row(i);
			Eigen::Vector3d URtmp = URi.row(i);
			Eigen::Vector3d ULtmp1 = ULi.row(i+1);
			Eigen::Vector3d URtmp1 = URi.row(i+1);

			Eigen::Vector3d ULbar = ULtmp1 + 0.5*(dt/dx)*(f(ULtmp1) - f(URtmp1)); //UL(i+1)
			Eigen::Vector3d URbar = URtmp + 0.5*(dt/dx)*(f(ULtmp) - f(URtmp));


			/*-------------------------------------------------------
			 * Solution of the piecewise constant Riemann problem pg 180
			 -------------------------------------------------------*/
			/*-------------------------------------------------------
			 * HLLC solver
			 -------------------------------------------------------*/
			Eigen::Vector3d hllcUL = URbar;
			Eigen::Vector3d hllcUR = ULbar;

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
				Pl = PressureScalar(hllcUL);
				Pr = PressureScalar(hllcUR);

				//velocity
				ul = hllcUL(1)/hllcUL(0);
				ur = hllcUR(1)/hllcUR(0);

				//soundspeed
				al = soundspeedScalar(hllcUL);
				ar = soundspeedScalar(hllcUR);



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
				//if (count == 0) std::cout << F.row(i) << std::endl;
		}
		//end of domain loop
		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > tstop) dt = tstop - t;
		t += dt;
		count += 1;

		if (count == 0) std::cout << dt << std::endl;
		//updating U
		for (int i=2; i<N+2; i++){
			//if (count == 0) std::cout << U.row(i) << '\t' << '\t';
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
			//if (count == 0) std::cout << U.row(i) << std::endl;
		}
		boundary_conditions(U);


	}while (t<tstop);
	std::cout << count << std::endl;
}

void output(MatrixN3 U, VectorN2 X, double tstop){

	std::stringstream filename;
	filename << "Euler_" << std::setw(3) << std::setfill('0') << tstop << "s.txt";

	std::ofstream outfile;
	outfile.open(filename.str().c_str());

	for (int i=0; i<(N+1); i++){

		double u = U(i+2, 1)/U(i+2, 0);
		double P = Pressure(U, i+2);
		//double e = (y-1)*(U(i+1, 2)-0.5*U(i+1, 1)*(U(i+1, 1)/U(i+1, 0)))/(U(i,0)*(y-1));
		double e = internalE(U, i+2);

		outfile << X(i) << '\t' << U(i+2, 0) << '\t' << u
				<< '\t' << P << '\t' << e << std::endl;
	}
	std::cout << "done" << std::endl;
}

int main(void){

	VectorN2 X; //mesh
	MatrixN3 U; //U(i)(0) = rho | U(i)(1) = rhou | U(i)(2) = E

	initial_conditions(X, U);
	MHM(U, tstop, CFL);
	output(U, X, 0.2);
}


