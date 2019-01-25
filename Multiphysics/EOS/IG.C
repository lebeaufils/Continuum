#include "IG.h"

double IdealGas::internalE(Eigen::MatrixXd U, int i){
	//E = ρ ( 0.5*V2 + e)
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double IdealGas::Pressure(Eigen::MatrixXd U, int i){
	double Pressure = (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0));
	return Pressure;
}

double IdealGas::PressureScalar(vector U){
	double Pressure = (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	return Pressure;
}

double IdealGas::soundspeed(Eigen::MatrixXd U, int i){
	double a = sqrt(y*(Pressure(U, i)/U(i, 0)));
	return a;
}

double IdealGas::soundspeedScalar(vector U){
	double a = sqrt(y*(PressureScalar(U)/U(0)));
	return a;
}

vector IdealGas::f(vector U){
	vector flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	flux(2) = (U(1)/U(0))*(U(2) + (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)));
	return flux;
}
