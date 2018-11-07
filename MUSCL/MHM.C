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
const double x0 = 0.3;
const double dx = L/N; //length of domain
const double CFL = 0.9;
const double y = 1.4; //heat capacity ratio Cp/Cv
const double tstop = 0.2;

const Eigen::Vector3d initialL(1.0, 0.75, 1.0); //density, velocity, pressure
const Eigen::Vector3d initialR(0.125, 0.0, 0.1); //left and right states
//const Eigen::Vector3d initialL(1.0, -2.0, 0.4);
//const Eigen::Vector3d initialR(1.0, 2.0, 0.4);
//const Eigen::Vector3d initialL(1.0, 0.0, 1000.0);
//const Eigen::Vector3d initialR(1.0, 0.0, 0.01);
//const Eigen::Vector3d initialL(5.99924, 19.5975, 460.894); //density, velocity, pressure
//const Eigen::Vector3d initialR(5.99242, -6.19633, 46.0950); //left and right states
//const Eigen::Vector3d initialL(1.0, -19.59745, 1000.0);
//const Eigen::Vector3d initialR(1.0, -19.59745, 0.01); // slowly moving contact discontinuities


typedef Eigen::Matrix<double, N+4, 3> MatrixN3;
typedef Eigen::Matrix<double, N+1, 3> MatrixN1;
typedef Eigen::Matrix<double, N+2, 1> VectorN2;
typedef Eigen::Matrix<double, N+1, 1> VectorN1;

void boundary_conditions(MatrixN3 &U){
	U.row(0) = U.row(1);
	U.row(1) = U.row(2); //muscl requires an extra ghost cell
	U.row(N) = U.row(N-1);
	U.row(N+1) = U.row(N);
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

	//E
	eulerL(2) = initialL(2)/(y-1) + 0.5*initialL(0)*initialL(1)*initialL(1);
	eulerR(2) = initialR(2)/(y-1) + 0.5*initialR(0)*initialR(1)*initialR(1);

	for (int i=0; i<N+1; i++){
		X(i) = i*dx;
		if (X(i)  < x0){
			U.row(i+2) = eulerL;
		}
		else U.row(i+2) = eulerR;
	}

	boundary_conditions(U);
}

/*
Eigen::Vector3d f(MatrixN2 U, int i){
	Eigen::Vector3d flux(U(i, 1),
			U(i, 1)*(U(i, 1)/U(i, 0))+(y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0)),
			(U(i, 1)/U(i, 0))*(U(i, 2) + (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0))));
	return flux;
}
*/

Eigen::Vector3d f(Eigen::Vector3d U){
	Eigen::Vector3d flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	flux(2) = (U(1)/U(0))*(U(2) + (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)));
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

	MatrixN1 F1;

	double dt = 0.01;
	double t = 0.0;
	int count = 0;
	do{

		for (int i=1; i<N+1; i++){

			Eigen::Vector3d diMinus = U.row(i) - U.row(i-1);
			Eigen::Vector3d diPlus = U.row(i+1) - U.row(i);
			Eigen::Vector3d di = 0.5*diMinus + 0.5*diPlus;

			/*-----------------------------------------------
			 * Slope limiter -- SuperBee
			 ----------------------------------------------*/

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
					epsilonR(j) = 2*2/(1+ri(j)); //need to refine beta using additional riemann problems
					epsiloni(j) = fmin(fmin(ri(j), epsilonR(j)), 2.0);
				}
			}

			Eigen::Vector3d diBar(0, 0, 0);
			for (int j=0; j<3; j++){
				 diBar(j) = epsiloni(j)*di(j);
			}
			//std::cout << diBar.transpose() << std::endl;


			/*-----------------------------------------------
			 * Slope limiter -- Van Leer
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
			*/

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
		for (int i=1; i<N; i++){

			Eigen::Vector3d ULtmp = ULi.row(i);
			Eigen::Vector3d URtmp = URi.row(i);
			Eigen::Vector3d ULtmp1 = ULi.row(i+1);
			Eigen::Vector3d URtmp1 = URi.row(i+1);

			Eigen::Vector3d ULbar = ULtmp1 + 0.5*(dt/dx)*(f(ULtmp1) - f(URtmp1));
			Eigen::Vector3d URbar = URtmp + 0.5*(dt/dx)*(f(ULtmp) - f(URtmp));

			/*-------------------------------------------------------
			 * Solution of the piecewise constant Riemann problem pg 180
			 -------------------------------------------------------*/
			/*-------------------------------------------------------
			 * HLLC solver
			 -------------------------------------------------------*/
			Eigen::Vector3d hllcUL = URbar;
			Eigen::Vector3d hllcUR = ULbar;

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
				Pl = (y-1)*(hllcUL(2) - 0.5*hllcUL(0)*pow((hllcUL(1)/hllcUL(0)),2.0));
				Pr = (y-1)*(hllcUR(2) - 0.5*hllcUR(0)*pow((hllcUR(1)/hllcUR(0)),2.0));

				//velocity
				ul = hllcUL(1)/hllcUL(0);
				ur = hllcUR(1)/hllcUR(0);

				//soundspeed
				al = sqrt(y*(Pl/dl));
				ar = sqrt(y*(Pr/dr));

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

		//updating U
		for (int i=2; i<N; i++){
			U.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
		}
		boundary_conditions(U);


		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > tstop) dt = tstop - t;
		t += dt;
		count += 1;

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
		double P = (y-1)*(U(i+2, 2) - 0.5*U(i+2, 0)*pow((U(i+2, 1)/U(i+2, 0)),2.0));
		//double e = (y-1)*(U(i+1, 2)-0.5*U(i+1, 1)*(U(i+1, 1)/U(i+1, 0)))/(U(i,0)*(y-1));
		double e = P/(U(i+2, 0)*(y - 1));

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

