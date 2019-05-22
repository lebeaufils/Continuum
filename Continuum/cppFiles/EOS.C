#include "../headerFiles/EOS.h"

std::shared_ptr<StateFunctions> StateFunctions::create(EOStype type){

	switch(type){
		case(EOS_IG):
			return std::make_shared<IdealGas>();
			break;
		case(EOS_SG):
			return std::make_shared<StiffenedGas>();
			break;
	}

}

IdealGas::IdealGas() : StateFunctions() {
	C = Cmatrix::Zero();
}

IdealGas::IdealGas(double y) : StateFunctions(y) {
	C = Cmatrix::Zero();
}

StiffenedGas::StiffenedGas() : StateFunctions(), Pref(0){
	C = Cmatrix::Zero();
}

StiffenedGas::StiffenedGas(double y, double Pref) : StateFunctions(y), Pref(Pref){
	C = Cmatrix::Zero();
}

//Obtaining internal energy (e) from conserved variable total E
double StateFunctions::internalE(matrix U, int i){
	//E = ρ ( 0.5*V2 + e)
	//note: this internalE is actually ρe 
	double e = U(i, 2)/U(i, 0) - 0.5*pow(U(i, 1)/U(i, 0), 2.0);
	return e;
}

double StateFunctions::internalE(vector4 U){
	//E = ρ ( 0.5*V2 + e)
	//note: this internalE is actually ρe 
	double e = U(2)/U(0) - 0.5*(pow(U(1)/U(0), 2.0) + pow(U(3)/U(0), 2));
	return e;
}

double StateFunctions::soundspeed(matrix U, int i){
	double a = sqrt(y*(Pressure(U, i)/U(i, 0)));
	return a;
}

double StateFunctions::soundspeed(vector U){
	double a = sqrt(y*(Pressure(U)/U(0)));
	return a;
}

double StateFunctions::soundspeed(vector4 U){
	double a = sqrt(y*(Pressure(U)/U(0)));
	return a;
}

vector StateFunctions::fluxes(vector U){
	vector flux;
	flux(0) = U(1);
	flux(1) = U(1)*(U(1)/U(0)) + Pressure(U);
	flux(2) = (U(1)/U(0))*(U(2) + Pressure(U));
	return flux;
}

vector4 StateFunctions::fluxes(vector4 U){ 
	//Here, U is either Ux or Uy
	vector4 flux;
	flux(0) = U(1); //momentum in sweep direction
	flux(1) = U(1)*(U(1)/U(0)) + Pressure(U);
	flux(2) = (U(1)/U(0))*(U(2) + Pressure(U));
	flux(3) = U(3)*U(1)/U(0);
	return flux;
}

/*----------------------------------------------------------------------------------
	Ideal Gas EOS
----------------------------------------------------------------------------------*/
double IdealGas::Pressure(matrix U, int i){
	if (U(i, 0) != 0){
		double Pressure = (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0));
		return Pressure;
	}
	else {
		return 0;
	}
}

double IdealGas::Pressure(vector U){
	if (U(0) != 0){
		double Pressure = (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0));
		return Pressure;
	}
	else {
		return 0;
	}
}

double IdealGas::Pressure(vector4 U){
	if (U(0) != 0){
		double Pressure = (y-1)*(U(2) - 0.5*U(0)*(pow(U(1)/U(0),2.0) + pow(U(3)/U(0), 2)));
		return Pressure;
	}
	else {
		return 0;
	}
}

/*double IdealGas::soundspeed(matrix U, int i){
	double a = sqrt(y*(Pressure(U, i)/U(i, 0)));
	return a;
}

double IdealGas::soundspeed(vector U){
	double a = sqrt(y*(Pressure(U)/U(0)));
	return a;
}*/

vector IdealGas::conservedVar(vector W){
	vector consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(1); //Density * Velocity
	consV(2) = W(2)/(y-1) + 0.5*W(0)*W(1)*W(1); //Total Energy
	return consV;
}

//if error, might be an incorrect arrangement of variables.
vector4 IdealGas::conservedVar2Dx(vector4 W){
	vector4 consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(1); //Density * Velocity
	consV(2) = W(3)/(y-1) + 0.5*W(0)*(pow(W(1), 2) + pow(W(2), 2)); //Total Energy
	consV(3) = W(0)*W(2); //Density * Velocity_y
	return consV;
}

vector4 IdealGas::conservedVar2Dy(vector4 W){
	vector4 consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(2); //Density * Velocity
	consV(2) = W(3)/(y-1) + 0.5*W(0)*(pow(W(1), 2) + pow(W(2), 2)); //Total Energy
	consV(3) = W(0)*W(1); //Density * Velocity_y
	return consV;
}

vector4 IdealGas::primitiveVar(vector4 U){
	//U is assumed to have the structure of conservedVar2Dx 
	vector4 primV;
	primV(0) = U(0);
	primV(1) = fmax(U(1)/U(0), 0);
	primV(2) = fmax(U(3)/U(0), 0);
	primV(3) = fmax(Pressure(U), 0);
	return primV;
}

/*----------------------------------------------------------------------------------
	Stiffened Gas EOS
----------------------------------------------------------------------------------*/

double StiffenedGas::Pressure(Eigen::MatrixXd U, int i){
	double Pressure = (y-1)*(U(i, 2) - 0.5*U(i, 0)*pow((U(i, 1)/U(i, 0)),2.0)) - y*Pref;
	return Pressure;
}

double StiffenedGas::Pressure(vector U){
	double Pressure = (y-1)*(U(2) - 0.5*U(0)*pow((U(1)/U(0)),2.0)) - y*Pref;
	return Pressure;
}

double StiffenedGas::Pressure(vector4 U){
	double Pressure = (y-1)*(U(2) - 0.5*U(0)*(pow(U(1)/U(0),2.0) + pow(U(3)/U(0), 2))) - y*Pref;
	return Pressure;
}

/*double StiffenedGas::soundspeed(Eigen::MatrixXd U, int i){
	double a = sqrt(y*((Pressure(U, i) + Pref)/U(i, 0)));
	return a;
}

double StiffenedGas::soundspeed(vector U){
	double a = sqrt(y*((Pressure(U) + Pref)/U(0)));
	return a;
}*/

vector StiffenedGas::conservedVar(vector W){
	vector consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(1); //Density * Velocity
	consV(2) = ((W(2)+y*Pref)/(y-1)) + 0.5*W(0)*W(1)*W(1); //Total Energy
	return consV;
}

//if error, might be an incorrect arrangement of variables.
vector4 StiffenedGas::conservedVar2Dx(vector4 W){
	vector4 consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(1); //Density * Velocity
	consV(2) = (W(3) + y*Pref)/(y-1) + 0.5*W(0)*(pow(W(1), 2) + pow(W(2), 2)); //Total Energy
	consV(3) = W(0)*W(2); //Density * Velocity_y
	return consV;
}

vector4 StiffenedGas::conservedVar2Dy(vector4 W){
	vector4 consV;
	consV(0) = W(0); //Density
	consV(1) = W(0)*W(2); //Density * Velocity
	consV(2) = (W(3) + y*Pref)/(y-1) + 0.5*W(0)*(pow(W(1), 2) + pow(W(2), 2)); //Total Energy
	consV(3) = W(0)*W(1); //Density * Velocity_y
	return consV;
}

vector4 StiffenedGas::primitiveVar(vector4 U){
	vector4 primV;
	primV(0) = U(0);
	primV(1) = U(1)/U(0);
	primV(2) = U(3)/U(0);
	primV(3) = Pressure(U);
	return primV;
}

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
	double Qk = sqrt(Ak/(P + Bk)); // Reciprocal of mass flux

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
	C(9) = W(2)*((y-1)/(y+1)) + Pref*(2*y)/(y+1);  //Bk (W(2)+y*Pref)*((y-1)/(y+1));
	C(10) = y;
	C(11) = W(0);
	C(12) = W(1);
	C(13) = W(2);
}

double StiffenedGas::fk(double pstar){
	//data-dependent constants
	double Ak = C(8);
	double Bk = C(9);
	double Qk = sqrt(Ak/(pstar + Bk));
	//double dstar = C(11)*((2*y*Pref + (y-1)*C(13) + (y+1)*pstar)/(2*y*Pref + (y+1)*C(13) + (y-1)*pstar));
	//double Qk = sqrt((C(13) - pstar)/(1./dstar - 1./C(11)));

	double Fk;

	if (pstar > C(13)){ //Shock
		Fk = (pstar - C(13))*Qk;
	}

	else { //Rarefraction
		Fk = C(4)*C(0)*(pow((pstar + Pref)/(C(13) + Pref), C(1)) - 1);
	}
	return Fk;
}

double StiffenedGas::fprimek(double pstar){
	double Ak = C(8);
	double Bk = C(9);
	double Qk = sqrt(Ak/(pstar + Bk)); 
	//double dstar = C(11)*((2*y*Pref + (y-1)*C(13) + (y+1)*pstar)/(2*y*Pref + (y+1)*C(13) + (y-1)*pstar));
	//double Qk = sqrt((C(13) - pstar)/(1./dstar - 1./C(11)));
	double Fkprime;

	if (pstar > C(13)){ //Shock
		Fkprime = Qk*(1 - (pstar - C(13))/(2.*(pstar + Bk))); 
	}

	else { //Rarefraction
		Fkprime = (1./(C(0)*C(11)))*pow((pstar + Pref)/(C(13) + Pref), -C(2));
	}
	return Fkprime;
}




