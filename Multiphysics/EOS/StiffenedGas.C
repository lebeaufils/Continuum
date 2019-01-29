#include "StiffenedGas.h"

void StiffenedGas::GetGamma(){
	std::cout << "Value of the Isentropic Expansion Factor, y" << std::endl;
	std::cin >> y;
	std::cin.ignore(32767, '\n');

	if (std::cin.fail()){
		std::cin.clear();
		std::cin.ignore(32767, '\n');
		std::cout << "Error, please enter a real number" << std::endl;
		GetGamma();
	}

	std::cout << "Value of the stiffening parameter P∞" << std::endl;
	std::cin >> y;
	std::cin.ignore(32767, '\n');

	if (std::cin.fail()){
		std::cin.clear();
		std::cin.ignore(32767, '\n');
		std::cout << "Error, please enter a real number" << std::endl;
		GetGamma();
	}
}

//Obtaining internal energy (e) from conserved variable total E
double StiffenedGas::internalE(Eigen::MatrixXd U, int i){
	//E = ρ ( 0.5*V2 + e)
	//note: this internalE is actually ρe 
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double StiffenedGas::Pressure(Eigen::MatrixXd U, int i){
	double Pressure = (y-1)*internalE(U, i) + y*Pref;
	return Pressure;
}

double StiffenedGas::PressureScalar(vector U){
	double Pressure = (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)) + y*Pref;
	return Pressure;
}

double StiffenedGas::soundspeed(Eigen::MatrixXd U, int i){
	double a = sqrt(y*(Pressure(U, i)/U(i, 0)));
	return a;
}

double StiffenedGas::soundspeedScalar(vector U){
	double a = sqrt(y*(PressureScalar(U)/U(0)));
	return a;
}

vector StiffenedGas::f(vector U){
	vector flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	flux(2) = (U(1)/U(0))*(U(2) + (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)));
	return flux;
}