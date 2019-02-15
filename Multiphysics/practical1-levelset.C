#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

const int N = 100; // number of cells
const double L = 1.0;
const double dx = L/N;
const double tstop = 0.2;

typedef Eigen::Matrix<double, N+2, 1> VectorN2;
typedef Eigen::Matrix<double, N, 1> VectorN;

void boundary_conditions(VectorN2 &U){
	U(0) = U(1);
	U(N+1) = U(N);
}

//first order upwind hamilton-jacobi method
void HJ_1(VectorN2 &phi, double velocity){

	double dt = (dx/abs(velocity));
	VectorN2 phi_1;

	double t = 0;
	int count = 0;
	do{

		for (int i=1; i<N+1; i++){

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
		boundary_conditions(phi);

		t += dt;
		std::cout << t << '\t' << dt << std::endl;
		count += 1;

	}while (t < tstop);
	std::cout << count << std::endl;

}

void forloopHJ(VectorN2 &phi, VectorN2 &phi_1, double velocity, double dt, int i){

	if (velocity <= 0){
		double tmp = phi(i) - (dt/dx)*velocity*(phi(i+1) - phi(i));
		phi_1(i) = tmp;
	}

	else if (velocity > 0){
		double tmp = phi(i) - (dt/dx)*velocity*(phi(i) - phi(i-1));
		phi_1(i) = tmp;
	}
}

void HJ_2(VectorN2 &phi, double velocity1, double velocity2){

	double dt1 = (dx/abs(velocity1));
	double dt2 = (dx/abs(velocity2));
	double dt = std::min(dt1, dt2);

	VectorN2 phi_1;

	VectorN X;
	for (int i=0; i<N; i++){
		X(i) = i*dx;
	}

	double t = 0;
	int count = 0;
	do{

		for (int i=1; i<N+1; i++){
			if (X(i-1) >= 0.5) forloopHJ(phi, phi_1, velocity1, dt, i);
			else forloopHJ(phi, phi_1, velocity2, dt, i);
		}
		phi = phi_1;
		boundary_conditions(phi);

		t += std::min(dt1, dt2);
		count += 1;

	}while (t < tstop);
	std::cout << count << std::endl;

}

int main(void){

	VectorN2 phi;
	VectorN X;

	//initial conditions
	for (int i=0; i<N; i++){
		X(i) = i*dx;
		phi(i+1) = X(i) - 0.5;
	}
	boundary_conditions(phi);


	//initial conditions 2(c)
	/*
	for (int i=0; i<N; i++){
		X(i) = i*dx;
		if (X(i) < 0.5){
			phi(i+1) = X(i) - 0.35;
			//catching underflow
			//if (abs(phi(i+1)) <= 1e-10) phi(i+1) = 0;
		}
		else phi(i+1) = -X(i) + 0.65;
	}
	boundary_conditions(phi);
	*/

	//initial output
	std::ofstream outfile_ini;
	outfile_ini.open("initial.txt");

	for (int i=0; i<N; i++){
		outfile_ini << X(i) << '\t' << phi(i+1) << std::endl;
		//std::cout << X(i) << '\t' << phi(i+1) << std::endl;
	}

	outfile_ini.close();

	//2(a)
	HJ_1(phi, 1.0);

	//2(b)
	//HJ_1(phi, -1.0);

	//2(c)
	//HJ_2(phi, 1.0, -1.0);


	//output
	std::ofstream outfile;
	outfile.open("data_levelset.txt");

	for (int i=0; i<N; i++){
		outfile << X(i) << '\t' << phi(i+1) << std::endl;
		//std::cout << X(i) << '\t' << phi(i+1) << std::endl;
	}

	outfile.close();

}

