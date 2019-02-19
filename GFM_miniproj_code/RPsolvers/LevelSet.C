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
	//phi(i+1) = X(i+1) - x0;
	double dist = abs(X(i+1) - x0);
	if (X(i+1) < x0) {
		phi(i+1) = -dist;
	}
	else {
		phi(i+1) = dist;
	}
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

void LevelSetFunction::signed_distance_function_1D_3(int i){
	// positive inside
	if (x0 > x1 || x1 > x2) {
		throw "Interface location incorrectly defined";
	}

	double dist = fmin(fmin(abs(X(i+1) - x0), abs(X(i+1) - x1)), abs(X(i+1) - x2));
	if (X(i+1) <= x0) {
		phi(i+1) = -dist;
	}

	else if (X(i+1) > x0 && X(i+1) <= x1){
		phi(i+1) = dist;
	}

	else if (X(i+1) > x1 && X(i+1) <= x2){
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

void LevelSetFunction::reinitialisation(){
//Sweeping in the positive x-direction
	for (int i=1; i<N+1; i++){
		//Consider the region phi > 0;
		if(phi(i) > 0){
			//Cell phi(i) is avaliable tp update if it is not adjacent to the interface
			if (phi(i+1) > 0 && phi(i-1) > 0){
				//Neither of its neighbours are <= 0
				//Taking the value of phi_x closest to the boudnary
				double phi_x = fmin(phi(i+1), phi(i-1));
				//if (i==N) phi_x = phi(i-1);
				//updating value based on the eikonal equation
				double phi_tmp = phi_x + dx; //taking the positive root
				//if (phi_tmp < phi(i)) phi(i) = phi_tmp;
				phi(i) = phi_tmp;
			}
		}

		else if (phi(i) < 0){
			//Consider the region phi < 0;
			if (phi(i+1) < 0 && phi(i-1) < 0){
				//Since phi is negative, the value closest to the boundary is the larger number
				double phi_x = fmax(phi(i+1), phi(i-1)); //taking the negative root
				//if (i==N) phi_x = phi(i-1);
				double phi_tmp = phi_x - dx;
				//if (phi_tmp > phi(i)) phi(i) = phi_tmp;
				phi(i) = phi_tmp;
			}
		}
	}
	boundary_conditions();
//Sweeping in the negative x-direction
	for (int i=N; i>0; i--){
		//Consider the region phi > 0;
		if(phi(i) > 0){
			//Cell phi(i) is avaliable tp update if it is not adjacent to the interface
			if (phi(i+1) > 0 && phi(i-1) > 0){
				//Neither of its neighbours are <= 0
				//Taking the value of phi_x closest to the boudnary
				double phi_x = fmin(phi(i+1), phi(i-1));
				//if (i==N) phi_x = phi(i-1);
				//updating value based on the eikonal equation
				double phi_tmp = phi_x + dx; //taking the positive root
				//if (phi_tmp < phi(i)) phi(i) = phi_tmp;
				phi(i) = phi_tmp;
			}
		}

		else if (phi(i) < 0){
			//Consider the region phi < 0;
			if (phi(i+1) < 0 && phi(i-1) < 0){
				//Since phi is negative, the value closest to the boundary is the larger number
				double phi_x = fmax(phi(i+1), phi(i-1)); //taking the negative root
				//if (i==N) phi_x = phi(i-1);
				double phi_tmp = phi_x - dx;
				//if (phi_tmp > phi(i)) phi(i) = phi_tmp;
				phi(i) = phi_tmp;
			}
		}
	}
	boundary_conditions();
}

/*void LevelSetFunction::reconstruction(){
	double newx0; 
	double newx1;

	int count = 0;
	for (int i=1; i<N+1; i++){
		int testsgn = (get_sgn(phi(i)) + get_sgn(phi(i+1)));
		if (count == 0){
			if (testsgn == 0) {
				double diff = dx*phi(i)/(phi(i+1) - phi(i));
				newx0 = phi(i) + diff;
				count = 1;
				i += 1;
			}
			else if (phi(i) == 0) {
				newx0 = phi(i);
				count = 1;
				i += 1;
			}
		}

		else if (count == 1) {
			if (testsgn == 0) {
				double diff = dx*phi(i)/(phi(i+1) - phi(i));
				newx1 = phi(i) + diff;
				count = 2;
			}

			else if (phi(i) == 0) {
				newx1 = phi(i);
				count = 2;
			}
		}
	}


	//std::cout << "Interface: " << newx0 << '\t' << newx1 << std::endl;
	// reconstructing the levelset function
	for (int i=0; i<N; i++){
		double dist = fmin(abs(X(i+1) - newx0), abs(X(i+1) - newx1));
		if (X(i+1) < newx0 || X(i+1) > newx1) {
			phi(i+1) = -dist;
		}

		else {
			phi(i+1) = dist;
		}
	}
}*/





