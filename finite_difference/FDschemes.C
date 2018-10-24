#include "FDschemes.H"

#include <fstream>
#include <iomanip>

//constructor
FiniteDifference::FiniteDifference(int n, double l, double a, double c)
	: N(n), L(l), a(a), c(c), X(N+2, 0), u(N+2, 0), u_1(N+2, 0), dx(l/n), dt(c*l/(n*a)){
}

LxF::LxF(int n, double l, double a, double c)
	: FiniteDifference(n, l, a, c){
}

LxW::LxW(int n, double l, double a, double c)
	: FiniteDifference(n, l, a, c){
}

WB::WB(int n, double l, double a, double c)
	: FiniteDifference(n, l, a, c){
}

void FiniteDifference::initial_conditions_step()
{
	for (int i=0; i<N; i++){
		X[i] = i*dx;

		if (X[i] < (L/2)){u[i+1] = 1.0;} //u starts from 1 to N+1 (index N)
		else {u[i] = 0.0;}
	}
	boundary_conditions();
}

void FiniteDifference::initial_conditions_square()
{
	for (int i=0; i<(N+1); i++){ //fixed boundary. check finite D
		X[i] = i*dx;

		if (X[i] < 0.3){u[i+1] = 0.0;} //u starts from 1 to N+1 (index N)
		else if (X[i] > 0.7){u[i+1] = 0.0;}
		else {u[i+1] = 1.0;}
	}
	boundary_conditions();
}

void FiniteDifference::boundary_conditions()
{
	u[0] = u[1];
	u[N+1] = u[N+2];
}

void FiniteDifference::output(std::string filename)
{
	std::ofstream outfile;
	outfile.open(filename);

	for (int i=1; i<(N+2); i++){
		outfile << X[i] << '\t' << u[i] << std::endl;
	}
}

void FiniteDifference::plotname(double t)
{
	std::stringstream filename;
	filename << "Plot_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}

//LxF
void LxF::lax_friedrichs()
{
	for (int i=1; i<(N+2); i++){
		u_1[i] = 0.5*(1+c)*u[i-1] + 0.5*(1-c)*u[i+1];
	}
	//updating the previous time-step
	for (int i=1; i<(N+2); i++){
		u[i] = u_1[i];
	}
}

void LxF::plotname(double t)
{
	std::stringstream filename;
	filename << "LxF_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}

//LxW
void LxW::lax_wendroff()
{
	for (int i=1; i<(N+2); i++){
		u_1[i] = 0.5*c*(1+c)*u[i-1] + (1-c*c)*u[i] - 0.5*c*(1-c)*u[i+1];
	}

	for (int i=1; i<(N+2); i++){
		u[i] = u_1[i];
	}
}

void LxW::plotname(double t)
{
	std::stringstream filename;
	filename << "LxW_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}

//WB
void WB::warming_beam()
{
	//note WB requires 2 ghost points to the left
	for (int i=2; i<(N+3); i++){
		u_1[i] = 0.5*c*(c-1)*u[i-2] + c*(2-c)*u[i-1] + 0.5*(c-1)*(c-2)*u[i];
	}

	for (int i=2; i<(N+3); i++){
			u[i] = u_1[i];
		}
}

void WB::initial_conditions_step()
{
	for (int i=0; i<(N+1); i++){
		X[i+2] = i*dx;

		if (X[i+2] < (L/2)){u[i+2] = 1.0;}
		else {u[i+2] = 0.0;}
	}
	boundary_conditions();
}

void WB::boundary_conditions()
{
	u[1] = u[2];
	u[0] = u[1];
}


void WB::output(std::string filename)
{
	std::ofstream outfile;
	outfile.open(filename);

	for (int i=2; i<(N+3); i++){
		outfile << X[i] << '\t' << u[i] << std::endl;
	}
}

void WB::plotname(double t)
{
	std::stringstream filename;
	filename << "WB_" << std::setw(3) << std::setfill('0') << t << "s.txt";
	output(filename.str().c_str());
}


