//ghost fluid part 1 practical

#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

const int N = 100; // number of cells
const double L = 4.0;
const double dx = L/N; //length of domain
const double x0 = 51*dx;
const double CFL = 0.9;
const double gamma1 = 1.4; //adiabatic exponent
const double gamma2 = 1.4; //material 2

//practial question 1
const Eigen::Vector3d initial1(2.0, 0.0, 9.8e5); //only change in density
const Eigen::Vector3d initial2(1.0, 0.0, 2.45e5);
//const Eigen::Vector3d initial1(1.0, 0.75, 1.0); //density, velocity, pressure
//const Eigen::Vector3d initial2(0.125, 0.0, 0.1); //left and right states
//onst Eigen::Vector3d initialL(1.0, -2.0, 0.4);
//const Eigen::Vector3d initialR(1.0, 2.0, 0.4);
//const Eigen::Vector3d initialL(1.0, 0.0, 1000.0);
//const Eigen::Vector3d initialR(1.0, 0.0, 0.01);
//const Eigen::Vector3d initialL(5.99924, 19.5975, 460.894); //density, velocity, pressure
//const Eigen::Vector3d initialR(5.99242, -6.19633, 46.0950); //left and right states
//const Eigen::Vector3d initialL(1.0, -19.59745, 1000.0);
//const Eigen::Vector3d initialR(1.0, -19.59745, 0.01); // slowly moving contact discontinuities


typedef Eigen::Matrix<double, N+2, 3> MatrixN2;
typedef Eigen::Matrix<double, N+1, 3> MatrixN1;
typedef Eigen::Matrix<double, N+2, 1> VectorN2;
typedef Eigen::Matrix<double, N+1, 1> VectorN1;

void boundary_conditions(MatrixN2 &U){
	//euler eqns
	U.row(0) = U.row(1);
	U.row(N+1) = U.row(N);
}

void boundary_levelset(VectorN2 &phi){
	//level set
	phi(0) = phi(1);
	phi(N+1) = phi(N);
}

void ghostBoundary(MatrixN2 Ureal, MatrixN2 &Ughost, int realmaterial, int i){ //where int realmaterial is the material index 1 or 2
	double yreal, yghost;
	double Preal, Pghost;
	double dghost;
	double dreal = Ureal(i, 0);
	double velocityreal = Ureal(i, 1)/Ureal(i, 0);
	double velocityghost;

	if (realmaterial == 1){
		yreal = gamma1;
		yghost = gamma2;
	}

	else {
		yreal = gamma2;
		yghost = gamma1;
	}

	Preal = (yreal-1)*(Ureal(i, 2) - 0.5*Ureal(i, 0)*pow((Ureal(i, 1)/Ureal(i, 0)),2.0)); //extrapolated Pressure
	Pghost = Preal;
	velocityghost = velocityreal;
	dghost = pow(pow(dreal, yreal)*(Pghost/Preal), 1.0/yghost);

	//std::cout << Pghost << '\t' << dghost << '\t' << velocityghost << std::endl;

	//setting the conservative form of the ghost variables
	Ughost(i, 0) = dghost;
	Ughost(i, 1) = velocityghost*dghost;
	Ughost(i, 2) = Pghost/(yghost-1) + 0.5*dghost*velocityghost*velocityghost;
}

void initial_conditions(VectorN2 &X, MatrixN2 &U1, MatrixN2 &U2, VectorN2 &phi){

	Eigen::Vector3d euler1;
	Eigen::Vector3d euler2;

	//rho
	euler1(0) = initial1(0);
	euler2(0) = initial2(0);

	//rhou
	euler1(1) = initial1(0)*initial2(1);
	euler2(1) = initial1(0)*initial2(1);

	//E
	euler1(2) = initial1(2)/(gamma1-1) + 0.5*initial1(0)*initial1(1)*initial1(1);
	euler2(2) = initial2(2)/(gamma2-1) + 0.5*initial2(0)*initial2(1)*initial2(1);


	for (int i=0; i<N+1; i++){
		X(i+1) = i*dx;
		if (X(i+1)  < x0){
			phi(i+1) = X(i+1) - x0; //signed distance func for 1D interface
			U1.row(i+1) = euler1; //setting real values
		}
		else{
			phi(i+1) = X(i+1) - x0;
			U2.row(i+1) = euler2; //setting real values
		}
	}

	//updating ghost values
	for (int i=0; i<N+1; i++){
		if (X(i+1)  < x0){
			ghostBoundary(U1, U2, 1, i+1); //setting ghost values through constant entropy extrapolation
		}
		else{
			ghostBoundary(U2, U1, 2, i+1);
		}
	}

	boundary_conditions(U1);
	boundary_conditions(U2);
	boundary_levelset(phi);
}

//sign function
int sgn(double x){
	return (x > 0) - (x < 0);
	/*
	Here's a more readable way to do it:

	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
	 */
}

void updateLevelSet(MatrixN2 &U1, MatrixN2 &U2, VectorN2 &phi, double dt){
	VectorN2 phi_1;
	double velocity;

	for (int i=1; i<N+1; i++){
		if (sgn(phi(i)) < 0){
			velocity = U1(i, 1)/U1(i, 0); //U1 is the 'left' material
		}

		else{
			velocity = U2(i, 1)/U2(i, 0); //U2 (right) material is the real material
		}

		if (velocity <= 0){
			double tmp = phi(i) - (dt/dx)*velocity*(phi(i+1) - phi(i));
			phi_1(i) = tmp;
		}

		else if (velocity > 0){
			double tmp = phi(i) - (dt/dx)*velocity*(phi(i) - phi(i-1));
			phi_1(i) = tmp;
		}
	}

	phi = phi_1;
	boundary_levelset(phi);
}

/*
double updateLevelSet(MatrixN2 &U, VectorN2 &phi, double dt){
	VectorN2 phi_1;

	for (int i=1; i<N+1; i++){
		double velocity = U(i, 1)/U(i, 0);
		if (velocity <= 0){
			double tmp = phi(i) - (dt/dx)*velocity*(phi(i+1) - phi(i));
			phi_1(i) = tmp;
		}

		else if (velocity > 0){
			double tmp = phi(i) - (dt/dx)*velocity*(phi(i) - phi(i-1));
			phi_1(i) = tmp;
		}
	}

	phi = phi_1;
	boundary_levelset(phi);

	int mid;
	for (int i=1; i<N+1; i++){
		//if the levelset function is at the zero contour
		if (sgn(phi(i)) == 0){
			mid = i;
		}

		//if the zero contour dosen't correspond to any discrete i value, take the midpoint
		else if (sgn(phi(i)) == -1 && sgn(phi(i+1)) == 1){
			mid = i;
		}

		else continue;
	}

	return (mid-1)*dx;
}
*/

/*
void reinitialise(VectorN2 &phi, VectorN2 X, double newx0){
	for (int i=0; i<N+1; i++){
		if (X(i)  < newx0){
			phi(i+1) = X(i) - newx0; //signed distance func for 1D interface
		}
		else{
			phi(i+1) = X(i) - newx0;
		}
	}
	boundary_levelset(phi);
}
*/

void computeFluxes(MatrixN2 U, MatrixN1 &F, double &Smax, const double y, int i){

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


	Eigen::Vector3d tmp;


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

	//soundspeed
	al = sqrt(y*(Pl/dl));
	ar = sqrt(y*(Pr/dr));

	/*
	 * pressure based wave speed estimate
	 */

	rhoavg = 0.5*(dl + dr);
	aavg = 0.5*(al + ar);

	Ppvrs = 0.5*(Pl + Pr) - 0.5*(ur - ul)*rhoavg*aavg;

	Pstar = fmax(0.0, Ppvrs);

	if (Pstar <= Pl){
		ql = 1.0;
	}

	else{
		ql = sqrt(1 + ((y+1)/(2*y))*((Pstar/Pl) - 1));
	}

	if (Pstar <= Pr){
		qr = 1.0;
	}

	else{
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

}

/*
void HLLCsolver(MatrixN2 &U, VectorN2 &phi, VectorN2 X, double tstop, double CFL){

	double oldx0 = x0;
	double Smax = 0;
	MatrixN1 F;

	//MatrixN2 U_1;

	double dt = 0.01;
	double t = 0.0;
	int count = 0;
	do{

		for (int i=0; i<N+1; i++){
			computeFluxes(U, F, Smax, i);
		}
		//end of domain loop

		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep

		//evolving level set function

		double newx0 = updateLevelSet(U, phi, dt);
		std::cout << newx0 << std::endl;
		if (newx0 > oldx0){
			reinitialise(phi, X, newx0);
			oldx0 = newx0;
		}

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
*/

void Ghostfluid_HLLCsolver(MatrixN2 &U1, MatrixN2 &U2, VectorN2 &phi, VectorN2 X, double tstop, double CFL){

	double Smax1 = 0;
	double Smax2 = 0;
	MatrixN1 F1;
	MatrixN1 F2;

	double dt = 0.0;
	double t = 0.0;
	int count = 0;
	do{
		//compute ghostfluid boundaries
		for (int i=0; i<N+1; i++){
			if (sgn(phi(i+1)) < 0) ghostBoundary(U1, U2, 1, i+1);
			else ghostBoundary(U1, U2, 1, i+1);
		}

		boundary_conditions(U1);
		boundary_conditions(U2);

		for (int i=0; i<N+1; i++){
			//looping across all real and ghost values
			computeFluxes(U1, F1, Smax1, gamma1, i);
			computeFluxes(U2, F2, Smax2, gamma2, i);
		}
		//find the maximum wave speed
		double Smax = std::max(Smax1, Smax2);
		//set timestep
		dt = CFL*(dx/Smax); //updates every timestep

		//evolving level set function
		//Real material velocity needs to be used.
		updateLevelSet(U1, U2, phi, dt);

		if (t + dt > tstop) dt = tstop - t;

		//updating U for material 1 and 2
		for (int i=1; i<N+1; i++){
			U1.row(i) = U1.row(i) - (dt/dx)*(F1.row(i) - F1.row(i-1));
			U2.row(i) = U2.row(i) - (dt/dx)*(F2.row(i) - F2.row(i-1));
		}
		boundary_conditions(U1);
		boundary_conditions(U2);


		t += dt;
		count += 1;
	}while(t < tstop);
	std::cout << count << std::endl;
}

void output(MatrixN2 U1, MatrixN2 U2, VectorN2 X, VectorN2 phi, double tstop){

	MatrixN2 U;
	double y;

	std::stringstream filename;
	filename << "Euler_" << std::setw(3) << std::setfill('0') << tstop << "s.txt";

	std::ofstream outfile;
	outfile.open(filename.str().c_str());

	for (int i=0; i<(N+1); i++){
		if (sgn(phi(i+1)) < 0){
			U = U1;
			y = gamma1;
		}

		else {
			U = U2;
			y = gamma2;
		}

		double u = U(i+1, 1)/U(i+1, 0);
		double P = (y-1)*(U(i+1, 2) - 0.5*U(i+1, 0)*pow((U(i+1, 1)/U(i+1, 0)),2.0));
		double e = P/(U(i+1, 0)*(y - 1));
		std::cout << u << '\t' << phi(i+1) << std::endl;

		outfile << X(i+1) << '\t' << U(i+1, 0) << '\t' << u
				<< '\t' << P << '\t' << e << '\t' << phi(i+1) << std::endl;
	}
	std::cout << "done" << std::endl;
}


int main(void){

	VectorN2 X; //mesh
	VectorN2 phi; //levelset function
	MatrixN2 U1; //U(i)(0) = rho | U(i)(1) = rhou | U(i)(2) = E
	MatrixN2 U2;

	initial_conditions(X, U1, U2, phi);
	double tstop = 0.0022;
	Ghostfluid_HLLCsolver(U1, U2, phi, X, tstop, CFL);
	output(U1, U2, X, phi, 0.2);
}







