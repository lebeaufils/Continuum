#include "../headerFiles/LevelSet.h"

int LevelSetMethods::get_sgn(double a){
	//sign function
	if (a > 0) return 1;
	else if (a < 0) return -1; 
	else return 0;
}


//----------------------------------------------------------------------------------------------------------
//1-Dimensional
//----------------------------------------------------------------------------------------------------------
void LevelSetMethods::initialise(LevelSet &ls, const Domain1D &domain){
	ls.phi.resize(domain.N+2, 1); //allows advection
}

void LevelSetMethods::boundary_conditions(LevelSet &ls, const Domain1D &domain){
	ls.phi(0) = ls.phi(1);
	ls.phi(domain.N+1) = ls.phi(domain.N);
}

void LevelSetMethods::signed_distance_function(LevelSet &ls, const Domain1D &domain, double x0){
	for (int i=0; i<domain.N; i++){
		ls.phi(i+1) = domain.X(i) - x0; //positive inside
	}
}

void LevelSetMethods::signed_distance_function(LevelSet &ls, const Domain1D &domain, double x0, double x1, int i){
	// positive inside
	if (x0 > x1) {
		throw "Interface location incorrectly defined";
	}

	double dist = fmin(abs(domain.X(i+1) - x0), abs(domain.X(i+1) - x1));
	if (domain.X(i+1) < x0 || domain.X(i+1) > x1) {
		ls.phi(i+1) = -dist;
	}

	else {
		ls.phi(i+1) = dist;
	}
}
//----------------------------------------------------------------------------------------------------------
//2-Dimensional
//----------------------------------------------------------------------------------------------------------
void LevelSetMethods::boundary_conditions(LevelSet &ls, const Domain2D &domain){
	//assigning ghost values in the x-direction 
	for (int j=0; j<domain.Ny; j++){
		ls.phi(0, j+1) = ls.phi(1, j+1);
		ls.phi(domain.Nx+3, j+1) = ls.phi(domain.Nx+2, j+1) = ls.phi(domain.Nx+1, j+1) = ls.phi(domain.Nx, j+1);
	} 
	//assigning ghost values in the y-direction
	for (int i=0; i<domain.Nx; i++){
		ls.phi(i+1, 0) = ls.phi(i+1, 1);
		ls.phi(i+1, domain.Ny+3) = ls.phi(i+1, domain.Ny+2) = ls.phi(i+1, domain.Ny+1) = ls.phi(i+1, domain.Ny);
	} 
}

void LevelSetMethods::initialise(LevelSet &ls, const Domain2D &domain, Polygon &poly){
	ls.phi = matrix::Zero(domain.Nx+4, domain.Ny+4);

	//Finding if the points are inside or outside the polygon
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (poly.point_in_polygon(domain.X(i, j))) ls.phi(i+1, j+1) = -1e6;
			else ls.phi(i+1, j+1) = 1e6;
		}
	}
	boundary_conditions(ls, domain);


	//The initial levelset has all points on the boundary set as 0
	for (int i=0; i<static_cast<int>(poly.surfacepoints.size()); i++){
		ls.phi(poly.surfacepoints[i].i+1, poly.surfacepoints[i].j+1) = 0;
	}

	fast_sweep(ls, domain);
}

void LevelSetMethods::initialise_circle(LevelSet &ls, Domain2D domain, double x0, double y0, double r){
	//This provides an exact levelset function
	ls.phi.resize(domain.Nx+4, domain.Ny+4);

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			//ls.phi(i, j) = pow(domain.X[i + j*domain.Nx][0] - x0, 2) + pow(domain.X[i + j*domain.Nx][1] - y0, 2) - pow(r, 2);
			ls.phi(i+1, j+1) = pow(domain.X(i, j).x - x0, 2) + pow(domain.X(i, j).y - y0, 2) - pow(r, 2);
		}
	}
	boundary_conditions(ls, domain);
	//ls.display_grid();

}

void LevelSetMethods::fast_sweep(LevelSet &ls, Domain2D domain){
	//Solver for the eikonal equation
	auto eikonal = [&ls, domain](int i, int j){
		//Consider the region phi > 0;
		if(ls.phi(i,j) > 0){
			//Cells that are not fixed (at the boundary) are avaliable for update
			//Taking the value of phi_x closest to the boudnary
			double phi_x = fmin(ls.phi(i+1, j), ls.phi(i-1, j));
			double phi_y = fmin(ls.phi(i, j+1), ls.phi(i, j-1));
			//updating value based on the eikonal equation
			double phi_tmp;
			if (pow(phi_x - phi_y, 2) >= pow(domain.dx, 2) + pow(domain.dy, 2)){
				phi_tmp = fmin(phi_x, phi_y) + sqrt(pow(domain.dx, 2) + pow(domain.dy, 2));
			}
			else {
				phi_tmp = (pow(domain.dy, 2)*phi_x + pow(domain.dx, 2)*phi_y + 
				domain.dx*domain.dy*sqrt(pow(domain.dx, 2) + pow(domain.dy, 2) - pow(phi_x - phi_y, 2)))/(pow(domain.dx, 2) + pow(domain.dy, 2));
			}
			if (phi_tmp < ls.phi(i,j)) ls.phi(i,j) = phi_tmp;
		}
		//Consider the region phi < 0;
		else if (ls.phi(i,j) < 0){
			//Since phi is negative, the value closest to the boundary is the larger number
			double phi_x = fmax(ls.phi(i+1, j), ls.phi(i-1, j));
			double phi_y = fmax(ls.phi(i, j+1), ls.phi(i, j-1));
			//if (i==N) phi_x = phi(i-1);
			double phi_tmp;
			if (pow(phi_x - phi_y, 2) >= pow(domain.dx, 2) + pow(domain.dy, 2)){
				phi_tmp = fmax(phi_x, phi_y) - sqrt(pow(domain.dx, 2) + pow(domain.dy, 2));
			}
			else {
				//taking the negative root
				phi_tmp = (pow(domain.dy, 2)*phi_x + pow(domain.dx, 2)*phi_y - 
				domain.dx*domain.dy*sqrt(pow(domain.dx, 2) + pow(domain.dy, 2) - pow(phi_x - phi_y, 2)))/(pow(domain.dx, 2) + pow(domain.dy, 2));
			}
			if (phi_tmp > ls.phi(i,j)) ls.phi(i,j) = phi_tmp;
		}
	};
	//Gauss Seidel iteration for alternate directions.
	//A total of four sweeps is performed over the entire computational domain
	//1) i = 1:I, j = 1:J
	//2) i = I:1, j = 1:J
	//3) i = I:1, J = J:1
	//4) i = 1:I, j = J:1

	//1) i = 1:I, j = 1:J
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			eikonal(i+1, j+1);
		}
	}

	//2) i = I:1, j = 1:J
	for (int i=domain.Nx-1; i>=0; i--){
		for (int j=0; j<domain.Ny; j++){
			eikonal(i+1, j+1);
		}
	}

	//3) i = I:1, J = J:1
	for (int i=domain.Nx-1; i>=0; i--){
		for (int j=domain.Ny-1; j>=0; j--){
			eikonal(i, j);
		}
	}

	//1) i = 1:I, j = 1:J
	for (int i=0; i<domain.Nx; i++){
		for (int j=domain.Ny-1; j>=0; j--){
			eikonal(i+1, j+1);
		}
	}

	for (int i=1; i<domain.Nx+1; i++){
		for (int j=1; j<domain.Ny+1; j++){
			std::cout << ls.phi(i,j) << '\t';
		}
		std::cout << std::endl;
	}
}

vector2 LevelSetMethods::normal(LevelSet ls, Domain2D domain, int i, int j){
	//Compute the normal vector using the central difference approximation
	double dphi_x = (ls.phi(i+1, j) - ls.phi(i-1, j))/(2*domain.dx);
	double dphi_y = (ls.phi(i, j+1) - ls.phi(i, j-1))/(2*domain.dy);
	double grad_phi = sqrt(dphi_x*dphi_x + dphi_y*dphi_y);

	vector2 n_i(dphi_x, dphi_y);
	if (grad_phi > 0){	
		n_i = n_i/grad_phi;
	}
	//else {
	//	throw "normal vector is 0";
	//}
	return n_i;
}

/*
void LevelSetFunction::signed_distance_function_1D(){
	//initial conditions
	for (int i=0; i<Test.N; i++){
		X(i+1) = i*dx;
		phi(i+1) = X(i) - x0;
	}
}*/

/*
//X is now not a member of the levelsetfucntion
void LevelSetFunction::signed_distance_function_1D(matrix X, int i, double x0){
	//phi(i+1) = X(i+1) - x0;
	double dist = abs(X(i+1) - x0);
	if (X(i+1) < x0) {
		phi(i+1) = -dist;
	}
	else {
		phi(i+1) = dist;
	}
}

void LevelSetFunction::signed_distance_function_1D_2(matrix X, int i, double x0, double x1){
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
	//comment this out
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
	//end
	// the above is for 3 material interfaces
	if (x1 > x2) {
		throw "Interface location incorrectly defined";
	}

	double dist = fmin(abs(X(i+1) - x1), abs(X(i+1) - x2));
	if (X(i+1) < x1 || X(i+1) > x2) {
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

void LevelSetFunction::reconstruction(){
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





