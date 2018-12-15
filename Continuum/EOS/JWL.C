#include "JWL.h"

double JWL::Pref(double d){ //function of density
	double dprime = d0/d;

	double Pref = A*exp(-R1*dprime) + B*exp(-R2*dprime);
	return Pref;
}

double JWL::eref(double d){
	double dprime = d0/d;

	double eref = A/(R1*d0)*exp(-R1*dprime) + B/(R2*d0)*exp(-R2*dprime);
	return eref;
}

double JWL::internalE(Eigen::MatrixXd U, int i){
	//E = ρ ( 0.5*V2 + e)
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double JWL::Pressure(Eigen::MatrixXd U, int i){
	double Pressure = Pref(U(i, 0)) + Gruneisen*U(i, 0)*(internalE(U, i) - eref(U(i, 0)));
	return Pressure;
}

double JWL::PressureScalar(vector U){
	double e = U(2)/U(0) - 0.5*pow(U(1)/U(0), 2.0);
	double Pressure = Pref(U(0)) + Gruneisen*U(0)*(e - eref(U(0)));
	return Pressure;
}

double JWL::soundspeed(Eigen::MatrixXd U, int i){
	double d = U(i, 0);
	double dprime = d0/d;
	double e = internalE(U, i);
	double csquare = (e - eref(d))*Gruneisen*(Gruneisen + 1)
			+ (1/pow(d, 2.0))*(A*R1*d0*pow(e, -R1*dprime) - B*R2*d0*pow(e, -R2*dprime));
	return sqrt(csquare);
}

double JWL::soundspeedScalar(vector U){
	double d = U(0);
	double dprime = d0/d;
	double e = U(2)/U(0) - 0.5*pow(U(1)/U(0), 2.0);
	double csquare = (e - eref(d))*Gruneisen*(Gruneisen + 1)
			+ (1/pow(d, 2.0))*(A*R1*d0*pow(e, -R1*dprime) - B*R2*d0*pow(e, -R2*dprime));
	return sqrt(csquare);
}
