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

	/*for(int i=0; i<domain.Nx; i++){
		for(int j=0; j<domain.Ny; j++){
			std::cout << ls.phi(i+1, j+1) << '\t';
		}
		std::cout << std::endl;
	}*/

	fast_sweep(ls, domain);

	/*for(int i=0; i<domain.Nx; i++){
		for(int j=0; j<domain.Ny; j++){
			std::cout << ls.phi(i+1, j+1) << '\t';
		}
		std::cout << std::endl;
	}*/

	std::ofstream outfile;
	outfile.open("extrapolation.txt");
	for(int i=0; i<domain.Nx; i++){
		for(int j=0; j<domain.Ny; j++){
			outfile << i*domain.dx << '\t' << j*domain.dy << '\t' << ls.phi(i+1, j+1) << std::endl;
		}
		outfile << std::endl;
	}
	outfile.close();
}

void LevelSetMethods::initialise_circle(LevelSet &ls, const Domain2D &domain, double x0, double y0, double r){
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

void LevelSetMethods::fast_sweep(LevelSet &ls, const Domain2D &domain){
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
				//if the quadratic solution is ill-defined
				phi_tmp = fmin(phi_x, phi_y) + sqrt(pow(domain.dx, 2) + pow(domain.dy, 2));
				//std::cout << "undefined " << i*domain.dx << '\t' << j*domain.dy << std::endl;
			}
			else {
				phi_tmp = (pow(domain.dy, 2)*phi_x + pow(domain.dx, 2)*phi_y + 
				domain.dx*domain.dy*sqrt(pow(domain.dx, 2) + pow(domain.dy, 2) - pow(phi_x - phi_y, 2)))/(pow(domain.dx, 2) + pow(domain.dy, 2));
				//std::cout << "          " << i*domain.dx << '\t' << j*domain.dy << std::endl;
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
				//std::cout << i*domain.dx << '\t' << j*domain.dy << std::endl;
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
			eikonal(i+1, j+1);
		}
	}

	//1) i = 1:I, j = 1:J
	for (int i=0; i<domain.Nx; i++){
		for (int j=domain.Ny-1; j>=0; j--){
			eikonal(i+1, j+1);
		}
	}

	/*for (int i=1; i<domain.Nx+1; i++){
		for (int j=1; j<domain.Ny+1; j++){
			std::cout << ls.phi(i,j) << '\t';
		}
		std::cout << std::endl;
	}*/

}

vector2 LevelSetMethods::normal(const LevelSet &ls, const Domain2D &domain, int i, int j){
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

double LevelSetMethods::interpolation_value(const LevelSet& ls, const Domain2D& domain, const Coordinates& p){
	//finding the 4 grid points surrounding point p
	int i = floor(p.x/domain.dx);
	int j = floor(p.y/domain.dy);

	//translation
	double x = p.x - domain.X(i, j).x;
	double y = p.y - domain.X(i, j).y;

	//bilinear interpolation
	double phi_xy = 0;
	for (int a=0; a<=1; a++){
		for (int b=0; b<=1; b++){
			//from the bilinear interpolation formnula,
			//if a is 0, contribution is (1-x), else if a is 1, contribution is x

			//additionally, to translate the interpolation square into the 
			phi_xy += ls.phi(a+i, b+j) * ((1-a)*(1-x) + a*x) * ((1-b)*(1-y) + b*y);
			//where phi(i, j) = c00, phi(i+1, j) = c10, phi(i, j+1) = c01 and phi(i+1, j+1) = c11
		}
	}
	return phi_xy;
}

vector2 LevelSetMethods::interpolation_gradient(const LevelSet& ls, const Domain2D& domain, const Coordinates& p){
	//differentiate the expression from interpolation value

	//As before,
	int i = floor(p.x/domain.dx);
	int j = floor(p.y/domain.dy);

	double x = p.x - domain.X(i, j).x;
	double y = p.y - domain.X(i, j).y;

	//using bilinear interpolation, taking the x and y derivatives
	double grad_x = 0;
	double grad_y = 0;
	for (int a=0; a<=1; a++){
		for (int b=0; b<=1; b++){
			grad_x += ls.phi(a+i, b+j) * (2*a-1) * ((1-b)*(1-y) + b*y);
			grad_y += ls.phi(a+i, b+j) * ((1-a)*(1-x) + a*x) * (2*b-1);
		}
	}
	return vector2(grad_x, grad_y);
}


//--------------------------------------------------------------------------------
// Calculating inertial properties using the levelsett function
//--------------------------------------------------------------------------------
double LevelSetMethods::smoothed_heaviside(const LevelSet &ls, const Domain2D &domain, int i, int j){
	const double e = 1.5*(domain.dx*domain.dy); //smoothness parameter
	const double pi = 3.14159265358979323846; //pi

	if (ls.phi(i, j) < -e){ //inside rigid body
		return 1;
	}
	else if (ls.phi(i, j) > e){ //outside rigid body
		return 0;
	}
	else { //within smoohed region
		return 0.5*(1 + ls.phi(i, j)/e + sin(pi*ls.phi(i, j)/e)/pi);
	}
}

double LevelSetMethods::mass(const LevelSet& ls, const Domain2D& domain, const Grain& gr){

	double m = 0;

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			m += smoothed_heaviside(ls,domain, i+1, j+1);
			//sums the number of cells within the rigid body, including the
			//transition zone
		}
	}
	m = gr.density*domain.dx*domain.dy*m; //rho * g^2 where g is the grid spacing

	return m;
}

vector2 LevelSetMethods::center_of_mass(const LevelSet& ls, const Domain2D& domain, const Grain& gr){
	//For a mass with density rho(r) within the solid, the integral of the weighted (density) position
	//coordinates relative to its center of mass is 0
	double c_x = 0;
	double c_y = 0;

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			c_x += smoothed_heaviside(ls,domain, i+1, j+1)*(domain.X(i, j).x);
			c_y += smoothed_heaviside(ls,domain, i+1, j+1)*(domain.X(i, j).y);
		}
	}
	double m = mass(ls, domain, gr);
	c_x = (gr.density*domain.dx*domain.dy/m)*c_x;
	c_y = (gr.density*domain.dx*domain.dy/m)*c_y;
	vector2 c(c_x, c_y);
	return c;
}

double LevelSetMethods::moment_of_inertia(const LevelSet& ls, const Domain2D& domain, const Grain& gr){
	//for rotation confined to a plane, the mass moment of inertia reduces to a scalar value.
	vector2 c = center_of_mass(ls, domain, gr);

	double inertia = 0;
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			inertia += smoothed_heaviside(ls,domain, i+1, j+1)*(i*domain.dx - c(0))*(j*domain.dy - c(1));
		}
	}
	inertia = -gr.density*domain.dx*domain.dy*inertia;
	return inertia;
}

//--------------------------------------------------------
// Motion
//--------------------------------------------------------
// Discrete equations of motion (assumes forces are provided (eg surface forces from fluid))

LevelSet LevelSetMethods::translation(const LevelSet& ls, const Domain2D& domain, const vector2& velocity){
	//x(n+1) = x(n) + v(t)*t
	//phi(x, t) = phi(x - vt, 0)
	//Temporary storage for levelset values at time t
	//The levelset can be updated using its initial values and current position
	LevelSet ls_t; //ls at time t
	ls_t.phi = matrix::Zero(domain.Nx+4, domain.Ny+4);

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates displacement(domain.X(i,j).x - velocity(0), domain.X(i, j).y - velocity(1));
			ls_t.phi(i+1, j+1) = interpolation_value(ls, domain, displacement);
		}
	}

	boundary_conditions(ls_t, domain);
	return ls_t;
}

void LevelSetMethods::rotation(LevelSet& ls, const Domain2D& domain){
	//the linear velocity is spatially dependent
	//vb = 2πω (r_rot × x)
}




