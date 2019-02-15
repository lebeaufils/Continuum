#include "EOS.h"

EOS::EOS()
	 : y(1.4){ //, Pref(0) {
		C = Eigen::Matrix<double, 14, 1>::Zero();
	 }

//Obtaining internal energy (e) from conserved variable total E
double EOS::internalE(Eigen::MatrixXd U, int i){
	//E = ρ ( 0.5*V2 + e)
	//note: this internalE is actually ρe 
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double EOS::soundspeed(Eigen::MatrixXd U, int i){
	double a = sqrt(y*(Pressure(U, i)/U(i, 0)));
	return a;
}

double EOS::soundspeedScalar(vector U){
	double a = sqrt(y*(PressureScalar(U)/U(0)));
	return a;
}

vector EOS::f(vector U){
	vector flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + PressureScalar(U);
	flux(2) = (U(1)/U(0))*(U(2) + PressureScalar(U));
	return flux;
}

void EOS::testing(){
	std::cout << C << std::endl;
}

 /*void EOS::assign_EOS_parameters(double gamma, double Pinfi){
 	y = gamma;
 	Pref = Pinfi;
 }*/

/*----------------------------------------------------------------------------------
	Ideal Gas EOS
----------------------------------------------------------------------------------*/

void IdealGas::GetGamma(){
	std::cout << "Value of the Isentropic Expansion Factor, y" << std::endl;
	std::cin >> y;
	std::cin.ignore(32767, '\n');

	if (std::cin.fail()){
		std::cin.clear();
		std::cin.ignore(32767, '\n');
		std::cout << "Error, please enter a real number" << std::endl;
		GetGamma();
	}
}

double IdealGas::Pressure(Eigen::MatrixXd U, int i){
	double Pressure = (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0));
	return Pressure;
}

double IdealGas::PressureScalar(vector U){
	double Pressure = (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
	return Pressure;
}

vector IdealGas::conservedVar(vector W){
	vector consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(1); //Density * Velocity
	consV(2) = W(2)/(y-1) + 0.5*W(0)*W(1)*W(1); //Total Energy
	return consV;
}

/*----------------------------------------------------------------------------------
	Stiffened Gas EOS
----------------------------------------------------------------------------------*/

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
	std::cin >> Pref;
	std::cin.ignore(32767, '\n');

	if (std::cin.fail()){
		std::cin.clear();
		std::cin.ignore(32767, '\n');
		std::cout << "Error, please enter a real number" << std::endl;
		GetGamma();
	}
}

double StiffenedGas::Pressure(Eigen::MatrixXd U, int i){
	double Pressure = (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0)) + y*Pref;
	return Pressure;
}

double StiffenedGas::PressureScalar(vector U){
	double Pressure = (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)) + y*Pref;
	return Pressure;
}

vector StiffenedGas::conservedVar(vector W){
	vector consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(1); //Density * Velocity
	consV(2) = (W(2)/(y-1) - y*Pref) + 0.5*W(0)*W(1)*W(1); //Total Energy
	return consV;
}

//-------------------------------------------------------------------------
//	Exact solver
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	Exact solver 
//-------------------------------------------------------------------------
void IdealGas::y_constants(vector W){
	C(0) = sqrt(y*W(2)/W(0));
	C(1) = (y-1)/(2*y); // y-1 / 2y
	C(2) = (y+1)/(2*y); // y+1 / 2y
	C(3) = (2*y)/(y-1);
	C(4) = 2./(y-1);
	C(5) = 2./(y+1);
	C(6) = (y-1)/(y+1);
	C(7) = (y-1)/2.;
	C(8) = (2./(y+1))/W(0);//Ak
	C(9) = W(2)*((y-1)/(y+1));  //Bk
	C(10) = y;
	C(11) = W(0);
	C(12) = W(1);
	C(13) = W(2);
}

double IdealGas::fk(double P){
	//data-dependent constants
	double Ak = C(8);
	double Bk = C(9);
	double Qk = sqrt(Ak/(P + Bk));

	//density = C(11);
	//velocity = C(12);
	//pressure = C(13);

	double Fk;

	if (P > C(13)){ //Shock
		Fk = (P - C(13))*Qk;
	}

	else { //Rarefraction
		Fk = C(4)*C(0)*(pow(P/C(13), C(1)) - 1);
	}
	return Fk;
}

double IdealGas::fprimek(double P){
	double Ak = C(8);
	double Bk = C(9);
	double Qk = sqrt(Ak/(P + Bk)); 
	double Fkprime;

	if (P > C(13)){ //Shock
		Fkprime = Qk*(1 - (P - C(13))/(2.*(P + Bk))); 
	}

	else { //Rarefraction
		Fkprime = (1./(C(11)*C(0)))*pow(P/C(13), -C(2));
	}
	return Fkprime;
}


void StiffenedGas::y_constants(vector W){
	C(0) = sqrt(y*(W(2) + Pref)/W(0));
	C(1) = (y-1)/(2*y); // y-1 / 2y
	C(2) = (y+1)/(2*y); // y+1 / 2y
	C(3) = (2*y)/(y-1);
	C(4) = 2./(y-1);
	C(5) = 2./(y+1);
	C(6) = (y-1)/(y+1);
	C(7) = (y-1)/2.;
	C(8) = (2./(y+1))/W(0);//Ak
	C(9) = W(2)*((y-1)/(y+1)) + Pref*(2*y)/(y+1);  //Bk
	C(10) = y;
	C(11) = W(0);
	C(12) = W(1);
	C(13) = W(2);
}

double StiffenedGas::fk(double P){
	//data-dependent constants
	double Ak = C(8);
	double Bk = C(9);
	double Qk = sqrt(Ak/(P + Bk));
	double Fk;

	if (P > C(13)){ //Shock
		Fk = (P - C(13))*Qk;
	}

	else { //Rarefraction
		Fk = C(4)*C(0)*(pow((P + Pref)/(C(13) + Pref), C(1)) - 1);
	}
	return Fk;
}

double StiffenedGas::fprimek(double P){
	double Ak = C(8);
	double Bk = C(9);
	double Qk = sqrt(Ak/(P + Bk)); 
	double Fkprime;

	if (P > C(13)){ //Shock
		Fkprime = Qk*(1 - (P - C(13))/(2.*(P + Bk))); 
	}

	else { //Rarefraction
		Fkprime = (1./(C(11)*C(0)))*pow((P + Pref)/(C(13) + Pref), -C(2));
	}
	return Fkprime;
}





