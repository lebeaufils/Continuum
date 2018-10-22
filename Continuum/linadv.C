#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

void initial_conditions(int N, double* X, double* u, double dx, double L){

	//Setting initial conditions
	for (int i=1; i<(N+1); i++){ //create domain from 0 to N+1, transmissive boundary
		X[i] = i*dx;

		if (X[i] < (L/2)){u[i] = 1.0;}
		else {u[i] = 0.0;}
	}
}

void boundary_conditions(int N, double* u){
	//copying boundary points to ghost points
	u[0] = u[1];
	u[N+1] = u[N];
}

void lax_friedrichs(int N, double* u, double* u_1, double c){
	//updating u[i], n+1
	for (int i=1; i<(N+2); i++){
		u_1[i] = 0.5*(1+c)*u[i-1] + 0.5*(1-c)*u[i+1];
	}
	//after looping over all the u points, u+1 becomes the previous timestep
	for (int i=1; i<(N+2); i++){
		u[i] = u_1[i];	
	}
}

void Lax_Wendroff(int N, double* u, double* u_1, double c){
	for (int i=1; i<(N+2); i++){
		u_1[i] = 0.5*c*(1+c)*u[i-1] + (1-c*c)*u[i] - 0.5*c*(1-c)*u[i+1];
	}
	for (int i=1; i<(N+2); i++){
		u[i] = u_1[i];
	}
}
/*
void Warming-Beam(int N, double* u, double* u_1, double c){
	for (int i=1; i<(N+2); i++){
		u_1[i] = 0.5*c*(c-1)*u[i-2] + 2*(2-c)*u[i-1] - 0.5*(c-1)*(c-2)*u[i];
	}
	for (int i=1; i<(N+2); i++){
		u[i] = u_1[i];
	}
}*/
// note the above scheme needs slightly different indexing

void output(string filename, int N, double *u, double *X){

	ofstream outfile;
	outfile.open(filename);

	for (int i=1; i<(N+2); i++){
			outfile << X[i] << '\t' << u[i] << endl;
	}
	outfile.close();
}

int main(void){
	//Model parameters
	int N = 100; //Number of cells
	double L = 1.0; //Domain length
	double a = 1.0; //Advection velocity
	double c = 0.9; //Courant (CFL) number

	//step size
	double dx = L/N;
	double dt = c*(dx/a);

	double *X = new double[N+2]; //creating domain
	double *u = new double[N+2]; //define state vector
	double *u_1 = new double[N+2]; //store previous timestep

	double t = 0;
	double tStop = 0.2;
	initial_conditions(N, X, u, dx, L);
	boundary_conditions(N, u);
	
	string initialfile = "initial_plot.txt";
	output(initialfile, N, u, X);

	int count = 0;
	do{
		if (t + dt > tStop){
			dt = tStop - t;
		}
		lax_friedrichs(N, u, u_1, c);
		boundary_conditions(N, u_1); //u_1 is the updated timestep u+1

		t += dt;
		count += 1;
		if (count%10 ==0){
			stringstream filename;
			filename << "LxF_" << setw(3) << setfill('0') << t << "s.txt";
			output(filename.str().c_str(), N, u_1, X);
			
		}

	}while (t < tStop);

	//output(N, u_1, X);
	cout << count << endl;

	delete[] X;
	delete[] u;
	delete[] u_1;
}
