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

void LevelSetMethods::initialise(LevelSet &ls, const Domain2D &domain, const Polygon &poly){
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
	ls.phi = matrix::Zero(domain.Nx+4, domain.Ny+4);

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
	if (p.x > domain.Lx || p.x < 0 || p.y > domain.Ly || p.y < 0){
		return 1e6; //returns a large positive value if domain is out of bounds
	}
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
			phi_xy += ls.phi(a+i+1, b+j+1) * ((1-a)*(1-x) + a*x) * ((1-b)*(1-y) + b*y);
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
// Utility functions
//--------------------------------------------------------------------------------
double LevelSetMethods::smoothed_heaviside(const LevelSet &ls, const Domain2D &domain, int i, int j){
	const double e = 1.5*domain.dx; //smoothness parameter
	const double pi = 3.14159265358979323846; //pi
	double phi_tmp = -ls.phi(i, j);

	if (phi_tmp < -e){ //inside rigid body
		return 0;
	}
	else if (phi_tmp > e){ //outside rigid body
		return 1;
	}
	else { //within smoohed region
		return 0.5*(1 + phi_tmp/e + sin(pi*phi_tmp/e)/pi);
	}
}

double LevelSetMethods::smoothed_delta(const LevelSet &ls, const Domain2D &domain, int i, int j){
	const double e = 1.5*domain.dx;//(domain.dx*domain.dy); //smoothness parameter
	const double pi = 3.14159265358979323846; //pi

	if (ls.phi(i, j) > e || ls.phi(i, j) < -e){ //inside rigid body
		return 0;
	}
	else { //outside rigid body
		return (1/(2*e) + 1/(2*e)*cos(pi*ls.phi(i, j)/e));
	}
}

LevelSet LevelSetMethods::merge(const std::vector<LevelSet>& levelsets, const Domain2D& domain){
	LevelSet ls(domain.Nx+4, domain.Ny+4);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			double min_phi = 1e6;
			for (int a=0; a<static_cast<int>(levelsets.size()); a++){
				if (levelsets[a].phi(i+1, j+1) < min_phi) min_phi = levelsets[a].phi(i+1, j+1);
			}
			ls.phi(i+1, j+1) = min_phi;
		}
	}
	boundary_conditions(ls, domain);

	return ls;
}

//LevelSet LevelSetMethods::intersection(const std::vector<LevelSet>&, const Domain2D& domain);
//LevelSet LevelSetMethods::difference(const std::vector<LevelSet>&, const Domain2D& domain);

//--------------------------------------------------------
// Forces
//--------------------------------------------------------
vector2 LevelSetMethods::force(const Euler2D& var, const LevelSet& ls, const Domain2D& domain){
//The force exerted on a solid (rigid body) by its surrounding liquid can be expressed as
//the integral sum of pressure normal to the surface area of the body
	vector2 f(0,0);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			vector2 n_i = LevelSetMethods::normal(ls, domain, i+1, j+1);
			double delta = LevelSetMethods::smoothed_delta(ls, domain, i+1, j+1);
			double p = var.state_function->Pressure(var.U(i+2, j+2)); //fluid pressure in cell i, j
			f += -(delta * p * n_i); //the normal direction is out of the rigid body
		}
	}
	return f;
}

double LevelSetMethods::torque (const Euler2D& var, const LevelSet& ls, const Domain2D& domain, const vector2& c){ 
//the total torque is the sum of cross products between the vector from the center of mass
	//to points lying on the interface (surface area) with the pressure normal to said area.
	//here, c is the position vector of the center of mass.
	double torque=0;
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			vector2 r_i(domain.X(i, j).x - c(0), domain.X(i, j).y - c(1));
			vector2 n_i = LevelSetMethods::normal(ls, domain, i+1, j+1);
			double delta = LevelSetMethods::smoothed_delta(ls, domain, i+1, j+1);
			double p = var.state_function->Pressure(var.U(i+2, j+2)); //fluid pressure in cell i, j
			
			//(x-c) x pn
			torque += -(delta * p * (r_i(0)*n_i(1) - r_i(1)*n_i(0)));//the normal direction is out of the rigid body
		}
	}
	return torque;
	//if this is combined with force into a subroutine, pressure can be computed only once
}

//--------------------------------------------------------
// Motion
//--------------------------------------------------------
// Discrete equations of motion (assumes forces are provided (eg surface forces from fluid))
// NOTE: in the paper by kawamoto, the levelset does not need to be updated as there is no evolution of surrounding fluid (GFM)
Coordinates LevelSetMethods::translation(const Coordinates& p, const vector2& v_c, double time){
	//x(n+1) = x(n) + v(t)*t
	//phi(x, t) = phi(x - vt, 0)

	Coordinates originalpos(p.x + v_c(0)*time, p.y + v_c(1)*time); //original position of the levelset
	return originalpos;
}

Coordinates LevelSetMethods::rotation(const Coordinates& p, const vector2& centroid, double f, double time){
	//ω is the angular velocity, which is a 'scalar' value in 2D,
	//f is the frequency of rotation, given by w/2pi
	//where a positive value corresponds to counter-clockwise rotation and negative clockwise.
	//the linear velocity is spatially dependent
	//vb = 2πω (r_rot × x)

	vector2 v = Rotor2::rotate_about(vector2(p.x, p.y), centroid, f, time); 
	//by reversing the rotation at point i, j, the "original position" of the level set is retrieved
	Coordinates originalpos(v(0), v(1)); //original position of the levelset
	return originalpos;
}

Coordinates LevelSetMethods::translation_reverse(const Coordinates& p, const vector2& v_c, double time){
	//x(n+1) = x(n) + v(t)*t
	//phi(x, t) = phi(x - vt, 0)

	Coordinates originalpos(p.x - v_c(0)*time, p.y - v_c(1)*time); //original position of the levelset
	return originalpos;
}

Coordinates LevelSetMethods::rotation_reverse(const Coordinates& p, const vector2& centroid, double f, double time){
	//ω is the angular velocity, which is a 'scalar' value in 2D,
	//f is the frequency of rotation, given by w/2pi
	//where a positive value corresponds to counter-clockwise rotation and negative clockwise.
	//the linear velocity is spatially dependent
	//vb = 2πω (r_rot × x)

	vector2 v = Rotor2::rotate_about(vector2(p.x, p.y), centroid, -f, time); 
	//by reversing the rotation at point i, j, the "original position" of the level set is retrieved
	Coordinates originalpos(v(0), v(1)); //original position of the levelset
	return originalpos;
}

LevelSet LevelSetMethods::motion(const LevelSet& ls, const Domain2D& domain, const vector2& centroid, const vector2& v_c, double w, double time){
	//Temporary storage for levelset values at time t
	//The levelset can be updated using its initial values and current position
	LevelSet ls_t; //ls at time t
	ls_t.phi = matrix::Zero(domain.Nx+4, domain.Ny+4);

	//simultaneous translation and rotation
	double pi = atan(1.0)*4;
	double frequency = w/(2*pi);

	for(int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates originalpos = domain.X(i, j);
			originalpos = translation_reverse(originalpos, v_c, time); //translate the point
			originalpos = rotation_reverse(originalpos, centroid, frequency, time); //rotate the point
			//double phi_tmp = interpolation_value(ls, domain, originalpos);
			ls_t.phi(i+1, j+1) = interpolation_value(ls, domain, originalpos);
			//if (phi_tmp < 1.5*domain.dx*domain.dy){ //within the rigidbody
			//	ls_t.phi(i+1, j+1) = phi_tmp;
			//}
			//else {
			//	ls_t.phi(i+1, j+1) = 1e6;
			//}
		}
	}

	fast_sweep(ls_t, domain);
	boundary_conditions(ls_t, domain);
	return ls_t;
}

//--------------------------------------------------------------
//Particles
//--------------------------------------------------------------
Particle::Particle(const Domain2D& domain, const Coordinates& center, double r) : ls(), centroid(0, 0), centre(0, 0), vc(0, 0), w(0), nodes(0), ref_nodes(0), miu(1), k_n(1e6), k_s(1e6), force(0, 0), torque(0) {
	LevelSetMethods::initialise_circle(ls, domain, center.x, center.y, r);
	centroid = center_of_mass(ls, *this, domain);
	centre = centroid;
	//std::cout << centroid.transpose() << std::endl;

	//seeding nodes for the circle
	//double circumference = 2*r*atan(1.0)*4;
	//double spacing = 2*r/10.;
	vector2 surface_p(0, 0);
	//find the first point on the zeroth levelset contour
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (ls.phi(i+1, j+1) == 0){
				surface_p(0) = domain.X(i, j).x;
				surface_p(1) = domain.X(i, j).y;
				break;
			}
		}
	}
	ref_nodes.push_back(surface_p);
	nodes.push_back(surface_p);
	//find the frequency of rotation
	int nodesize = 32; //circumference/spacing = 10pi
	double t = 0;
	for (int dt=1; dt<nodesize; dt++){
		vector2 v = Rotor2::rotate_about(surface_p, vector2(center.x, center.y), 2*3.14159/nodesize, t+dt); 
		ref_nodes.push_back(v);
		nodes.push_back(v);
	}			
}

Particle::Particle(const Polygon& poly, const Domain2D& domain) : ls(), centroid(0, 0), centre(0, 0), vc(0, 0), w(0), nodes(0), ref_nodes(0), miu(1), k_n(1000), k_s(1000), force(0, 0), torque(0) {
	LevelSetMethods::initialise(ls, domain, poly);
	centroid = center_of_mass(ls, *this, domain);
	centre = centroid;

	double circumference = static_cast<int>(poly.surfacepoints.size())*domain.dx;
	double diameter = circumference/3.2;
	int spacing = floor(diameter/(10*domain.dx));
	for (int a=0; a < static_cast<int>(poly.surfacepoints.size()); a+=spacing) {
		Coordinates tmp = domain.X(poly.surfacepoints[a].i, poly.surfacepoints[a].j);
		ref_nodes.push_back(vector2(tmp.x, tmp.y));
		nodes.push_back(vector2(tmp.x, tmp.y));
	}
	//for (int a=0; a < static_cast<int>(nodes.size()); a++) {
		//nodes[a].display();
	//}
}

Particle::Particle(const Particle& gr) : ls(gr.ls), centroid(gr.centroid), centre(gr.centre), vc(gr.vc), w(gr.w), nodes(0), ref_nodes(0), miu(gr.miu), k_n(gr.k_n), k_s(gr.k_s), force(gr.force), torque(gr.torque) {
	 for (int a=0; a<static_cast<int>(gr.nodes.size()); a++){
	 	nodes.push_back(gr.nodes[a]);
	 	ref_nodes.push_back(gr.ref_nodes[a]);
	 }
}
void Particle::set_velocity(const vector2& trans_velocity, double angular_velocity){
	vc = trans_velocity;
	w = angular_velocity;
}

double Particle::mass(const Particle& gr, const Domain2D& domain){

	double m = 0;

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			m += LevelSetMethods::smoothed_heaviside(gr.ls, domain, i+1, j+1);
			//sums the number of cells within the rigid body, including the
			//transition zone
		}
	}
	m = gr.density*domain.dx*domain.dy*m; //rho * g^2 where g is the grid spacing

	return m;
}

vector2 Particle::center_of_mass(const LevelSet& ls, const Particle& gr, const Domain2D& domain){
	//For a mass with density rho(r) within the solid, the integral of the weighted (density) position
	//coordinates relative to its center of mass is 0
	double c_x = 0;
	double c_y = 0;

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			c_x += LevelSetMethods::smoothed_heaviside(ls, domain, i+1, j+1)*(domain.X(i, j).x);
			c_y += LevelSetMethods::smoothed_heaviside(ls, domain, i+1, j+1)*(domain.X(i, j).y);
		}
	}
	double m = mass(gr, domain);
	c_x = (gr.density*domain.dx*domain.dy/m)*c_x;
	c_y = (gr.density*domain.dx*domain.dy/m)*c_y;
	vector2 c(c_x, c_y);
	return c;
}

double Particle::moment_of_inertia(const LevelSet& ls, const Particle& gr, const Domain2D& domain){
	//for rotation confined to a plane, the mass moment of inertia reduces to a scalar value.
	vector2 c = center_of_mass(ls, gr, domain);
	//in 2d, moment of inertia is a 1x1 tensor given by J = sum(m_i * (x_i^2 + y_i^2))
	double inertia = 0;
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			inertia += LevelSetMethods::smoothed_heaviside(ls, domain, i+1, j+1)*(pow((i*domain.dx - c(0)),2) + pow((j*domain.dy - c(1)),2));
		}
	}
	inertia = -gr.density*domain.dx*domain.dy*inertia;
	return inertia;
}

vector2 Particle::velocity(const Coordinates& p, const Particle& gr){
	vector2 v(0, 0);
	vector2 r(p.x-gr.centroid(0), p.y-gr.centroid(1));
	//v = vc + w cross r
	v(0) = gr.vc(0) - r(1)*gr.w;
	v(1) = gr.vc(1) + r(0)*gr.w;
	return v;
}

LevelSet Particle::merge(const std::vector<Particle>& particles, const Domain2D& domain){
	LevelSet ls(domain.Nx+4, domain.Ny+4);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			double min_phi = 1e6;
			for (int a=0; a<static_cast<int>(particles.size()); a++){
				if (particles[a].ls.phi(i+1, j+1) < min_phi) min_phi = particles[a].ls.phi(i+1, j+1);
			}
			ls.phi(i+1, j+1) = min_phi;
		}
	}
	LevelSetMethods::boundary_conditions(ls, domain);

	return ls;
}

LevelSet Particle::motion(const Domain2D& domain, double time){
	//Temporary storage for levelset values at time t
	//The levelset can be updated using its initial values and current position
	LevelSet ls_t; //ls at time t
	ls_t.phi = matrix::Zero(domain.Nx+4, domain.Ny+4);

	//simultaneous translation and rotation
	double pi = atan(1.0)*4;
	double frequency = w/(2*pi);
	//std::cout << "frequency = " << frequency << std::endl;

	for(int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates originalpos = domain.X(i, j);
			//originalpos.display();
			//originalpos = LevelSetMethods::rotation(originalpos, centroid, frequency, time); //rotate the point			
			originalpos = LevelSetMethods::rotation_reverse(originalpos, centroid, frequency, time); //rotate the point
			originalpos = LevelSetMethods::translation_reverse(originalpos, vc, time); //translate the point
			//originalpos.display();
			//originalpos.display();
			//double phi_tmp = interpolation_value(ls, domain, originalpos);
			ls_t.phi(i+1, j+1) = LevelSetMethods::interpolation_value(ls, domain, originalpos);
			//std::cout << "level set diff = " << ls_t.phi(i+1, j+1) << '\t' << ls.phi(i+1, j+1) << std::endl;
			//std::cout << std::endl;
		}
		//std::cout << std::endl;
	}
	//move the nodes from ref, store a seperate node list
	for (int a=0; a<static_cast<int>(nodes.size()); a++){
		//std::cout << nodes[a].transpose() << std::endl;
		Coordinates node_a(ref_nodes[a]); //always start from reference nodes
		node_a = LevelSetMethods::rotation(node_a, centroid, frequency, time);
		node_a = LevelSetMethods::translation(node_a, vc, time); 
		nodes[a] = vector2(node_a.x, node_a.y);
	}
////////
	LevelSetMethods::fast_sweep(ls_t, domain);
	LevelSetMethods::boundary_conditions(ls_t, domain);
	return ls_t;
}

vector2 Particle::cross(double w, const vector2& pos){
 	return vector2(-w*pos(1), w*pos(0));
}

double Particle::cross(const vector2& m, const vector2& n){
	return m(0)*n(1) - m(1)*n(0);
}

///////
void Moving_RB::add_particle(const Polygon& poly, const Domain2D& domain){
	//Particle new_particle(poly, domain);
	particles.push_back(Particle(poly, domain));
}

void Moving_RB::add_sphere(const Domain2D& domain, const Coordinates& center , double r){
	//Particle new_particle(domain, center, r);
	particles.push_back(Particle(domain, center, r));
}



