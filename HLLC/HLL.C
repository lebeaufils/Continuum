#include <iostream>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense> //dense matrices only. use sparse if needed

const int N = 100; // number of cells
const double L = 1.0;
const double x0 = 0.3;
const double dx = L/N; //length of domain
const double CFL = 0.9;
const double y = 1.4; //heat capacity ratio Cp/Cv
const Eigen::Vector3d initialL(1.0, 0.75, 1.0); //density, velocity, pressure
const Eigen::Vector3d initialR(0.125, 0.0, 0.1); //left and right states

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
	eulerL(2) = initialL(2)/(y-1) + 0.5*initialL(0)*initialL(1)*initialL(1);
	eulerR(2) = initialR(2)/(y-1) + 0.5*initialR(0)*initialR(1)*initialR(1);

	for (int i=0; i<N+2; i++){
		X(i) = i*dx;
		if (X(i)  < x0){
			U.row(i) = eulerL;
		}
		else U.row(i) = eulerR;
	}

	boundary_conditions(U);

}

//primitive variables
double pressure(MatrixN2 U, int i)
{
	//return (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow(((U(i, 1)/U(i, 0)),2.0)));
	//return (y-1)*(U(i,2)-0.5*U(i,1)/U(i,0)); //E-(0.5*(v^2))*rho*(y-1)
	return (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2));
}

double density(MatrixN2 U, int i)
{
	return U(i,0);
}

double momentum(MatrixN2 U, int i)
{
	return U(i,1);
}

double energy(MatrixN2 U, int i)
{
	return U(i,2);
}

double velocity(MatrixN2 U, int i)
{
	return U(i,1)/U(i,0);
}

double internalE(MatrixN2 U, int i)
{
	return (y-1)*(U(i,2)-0.5*U(i,1)*velocity(U, i))/(U(i,0)*(y-1));
}

//equation flux
double fluxmass(double density, double velocity)
{
	return density*velocity;
}

double fluxmom(double momentum, double velocity, double pressure)
{
	return momentum*velocity + pressure;
}

double fluxenergy(double energy, double velocity, double pressure)
{
	return velocity*(energy + pressure);
}

double computeSoundSpeed(double pressure, double density){
	return sqrt(y*(pressure/density));
}

double computeMaxWaveSpeed(double a, double velocity){
	return a + abs(velocity); //where a is sound speed
}

//calculate star state
void HLLsolver(MatrixN2 &U, double tstop, double CFL){

	double al;
	double ar;
	double ul; //note u is the eigenvalue of the euler equation
	double ur;
	double SL; //where S is the speed of acoustic waves x/t
	double SR;
	double Smax = 0;
	MatrixN1 F;
	Eigen::Vector3d FL;
	Eigen::Vector3d FR;
	Eigen::Vector3d tmp;

	MatrixN2 U_1;

	double dt = 0.01;
	double t = 0;
	int count = 0;
	do{

		for (int i=0; i<N+1; i++){
			al = computeSoundSpeed(pressure(U, i), density(U, i));
			ar = computeSoundSpeed(pressure(U, i+1), density(U, i+1));

			ul = velocity(U, i);
			ur = velocity(U, i+1);

			SR = std::max(al + abs(ul), ar + abs(ur));
			SL = -SR;

			FL << fluxmass(density(U, i), velocity(U, i)),
					fluxmom(momentum(U, i), velocity(U, i), pressure(U, i)),
					fluxenergy(energy(U, i), velocity(U, i), pressure(U, i));
			FR << fluxmass(density(U, i+1), velocity(U, i+1)),
							fluxmom(momentum(U, i+1), velocity(U, i+1), pressure(U, i+1)),
							fluxenergy(energy(U, i+1), velocity(U, i+1), pressure(U, i+1));
			//std::cout << FR.transpose() << std::endl;

			if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

			if (0 <= SL){
				F.row(i) = FL;
			}

			//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
			else if (SL <= 0 <= SR){
				//Rusanov Flux
				//F.row(i) = 0.5*(FL + FR) - 0.5*SR*(U.row(i+1) - U.row(i)); //matrices of diff size
				tmp = (U.row(i+1) - U.row(i));
				F.row(i) = (SR*FL - SL*FR + SL*SR*tmp)/(SR - SL);
			}

			else if (0 >= SR){
				F.row(i) = FR;
			}
		}

		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > tstop) dt = tstop - t;

		for (int i=1; i<N+1; i++){
			U_1.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
		}

		//updating U
		U = U_1;
		boundary_conditions(U);

		t += dt;
		count += 1;
		if (count%50 == 0) std::cout << dt << std::endl << std::endl;

	}while(t < tstop);

	std::cout << count << std::endl;
}

//calculate star state
void HLLpressuresolver(MatrixN2 &U, double tstop, double CFL){

	double al;
	double ar;
	double ul; //note u is the eigenvalue of the euler equation
	double ur;
	double pl; //pressure
	double pr;
	double dl; //density
	double dr;
	double SL; //where S is the speed of acoustic waves x/t
	double SR;
	double Smax = 0;
	double Ppvrs;
	double Pstar;
	double rhoavg;
	double aavg;
	double ql;
	double qr;
	MatrixN1 F;
	Eigen::Vector3d FL;
	Eigen::Vector3d FR;
	Eigen::Vector3d tmp;

	MatrixN2 U_1;

	double dt = 0.01;
	double t = 0;
	int count = 0;
	do{

		for (int i=0; i<N+1; i++){
			al = computeSoundSpeed(pressure(U, i), density(U, i));
			ar = computeSoundSpeed(pressure(U, i+1), density(U, i+1));

			ul = velocity(U, i);
			ur = velocity(U, i+1);
			pl = pressure(U, i);
			pr = pressure(U, i+1);
			dl = density(U, i);
			dr = density(U, i+1);

			rhoavg = 0.5*(dl + dr);
			aavg = 0.5*(al + ar);

			Ppvrs = 0.5*(pl + pr) - 0.5*(ur - ul)*rhoavg*aavg;

			//std::cout << rhoavg << '\t' << aavg << '\t' << Ppvrs << '\t' << dr << std::endl;

			//if (Ppvrs >= 0) Pstar = Ppvrs;
			//else if (Ppvrs < 0) Pstar = 0;
			Pstar = fmax(0.0, Ppvrs);

			if (Pstar <= pl){
				ql = 1.0;
			}

			else if (Pstar > pl){
				ql = sqrt(1 + ((y+1)/(2*y))*((Pstar/pl) - 1));
			}

			if (Pstar <= pr){
				qr = 1.0;
			}

			else if (Pstar > pr){
				qr = sqrt(1 + ((y+1)/(2*y))*((Pstar/pr) - 1));
			}

			SR = ur + ar*qr;
			SL = ul - al*ql;

			//std::cout << SR << '\t' << SL << std::endl;

			FL << fluxmass(density(U, i), velocity(U, i)),
					fluxmom(momentum(U, i), velocity(U, i), pressure(U, i)),
					fluxenergy(energy(U, i), velocity(U, i), pressure(U, i));
			FR << fluxmass(density(U, i+1), velocity(U, i+1)),
							fluxmom(momentum(U, i+1), velocity(U, i+1), pressure(U, i+1)),
							fluxenergy(energy(U, i+1), velocity(U, i+1), pressure(U, i+1));
			//std::cout << FR.transpose() << std::endl;

			if (std::max(abs(SR), abs(SL)) > Smax) Smax = std::max(abs(SR), abs(SL));

			if (0 <= SL){
				F.row(i) = FL;
			}

			//non trivial, subsonic case, region between the 2 acoustic waves. Fhll.
			else if (SL <= 0 <= SR){
				//Rusanov Flux
				//F.row(i) = 0.5*(FL + FR) - 0.5*SR*(U.row(i+1) - U.row(i)); //matrices of diff size
				tmp = (U.row(i+1) - U.row(i));
				F.row(i) = (SR*FL - SL*FR + SL*SR*tmp)/(SR - SL);
			}

			else if (0 >= SR){
				F.row(i) = FR;
			}
		}

		dt = CFL*(dx/Smax); //updates every timestep
		if (t + dt > tstop) dt = tstop - t;

		for (int i=1; i<N+1; i++){
			U_1.row(i) = U.row(i) - (dt/dx)*(F.row(i) - F.row(i-1));
		}

		//updating U
		U = U_1;
		boundary_conditions(U);

		t += dt;
		count += 1;
		if (count%50 == 0) std::cout << dt << std::endl << std::endl;

	}while(t < tstop);

	std::cout << count << std::endl;
}


void output(MatrixN2 U, VectorN2 X, double tstop){

	std::stringstream filename;
	filename << "Euler_" << std::setw(3) << std::setfill('0') << tstop << "s.txt";

	std::ofstream outfile;
	outfile.open(filename.str().c_str());

	for (int i=0; i<(N+1); i++){

		outfile << X(i) << '\t' << density(U, i+1) << '\t' << velocity(U, i+1)
				<< '\t' << pressure(U, i+1) << '\t' << internalE(U, i+1) << std::endl;
	}
	std::cout << "done" << std::endl;
}


int main(void){

	VectorN2 X; //mesh
	MatrixN2 U; //U(i)(0) = rho | U(i)(1) = rhou | U(i)(2) = E

	initial_conditions(X, U);

	double tstop = 0.1;

	//HLLsolver(U, tstop, CFL);
	HLLpressuresolver(U, tstop, CFL);
	//HLLCsolver(U, tstop, CFL);
	std::cout << "done" << std::endl;
	output(U, X, 0.2);
}

