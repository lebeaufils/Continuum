#include "LevelSet.h"

LevelSetFunction::LevelSetFunction(gfmTests Test)
	: N(Test.N), X(N+2, 1), dx(Test.L/Test.N), x0(Test.x0), x1(Test.x1), x2(Test.x2), phi(N+2, 1), sgn(0) {}

int LevelSetFunction::get_sgn(double a){
	//sign function
	if (a > 0) return 1;
	else if (a < 0) return -1; 
	else return 0;
}

void LevelSetFunction::boundary_conditions(){
	phi(0) = phi(1);
	phi(N+1) = phi(N);
}

/*
void LevelSetFunction::signed_distance_function_1D(){
	//initial conditions
	for (int i=0; i<Test.N; i++){
		X(i+1) = i*dx;
		phi(i+1) = X(i) - x0;
	}
}*/

void LevelSetFunction::signed_distance_function_1D(int i){
	phi(i+1) = X(i+1) - x0;
}

void LevelSetFunction::signed_distance_function_1D_2(int i){
	// positive inside
	if (x0 > x1) {
		throw "Interface location incorrectly defined";
	}

	double dist = fmin(abs(X(i+1) - x0), abs(X(i+1) - x1));
	if (X(i+1) < x0 || X(i+1) > x1) {
		phi(i+1) = -dist;
	}

	else {
		phi(i+1) = dist;
	}
}

double LevelSetFunction::HJ_FirstOrder(double velocity, double dt, int i){ //velocity and dt to follow numerical method

	double phi_1;

	if (velocity <= 0){
		phi_1 = phi(i) - velocity*(dt/dx)*(phi(i+1) - phi(i));
	}
	else {
		phi_1 = phi(i) - velocity*(dt/dx)*(phi(i) - phi(i-1));
	}

	return phi_1;
}

