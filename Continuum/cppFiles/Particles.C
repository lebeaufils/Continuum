#include "../headerFiles/Particles.h"
//--------------------------------------------------------------
//Particles
//--------------------------------------------------------------
Particle::Particle(const Domain2D& domain, const Coordinates& center, double r) : ls(), centroid(0, 0), centre(0, 0), size(0), vc(0, 0), w(0), s(0, 0), theta(0), nodes(0), ref_nodes(0), damping_coefficient(0.24), miu(0.26), k_n(1e5), k_s(1e5), force(0, 0), force_t(0, 0), torque(0), in_collision(false) {
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
			if (ls.phi(i+domain.buffer, j+domain.buffer) == 0){
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
	//calculate the size of the bounding box
	size = AABB_extent(ls);			
}

Particle::Particle(const Polygon& poly, const Domain2D& domain) : ls(), centroid(0, 0), centre(0, 0), size(0), vc(0, 0), w(0), s(0, 0), theta(0), nodes(0), ref_nodes(0), damping_coefficient(0.24), miu(0.26), k_n(1e5), k_s(1e5), force(0, 0), force_t(0, 0), torque(0), in_collision(false) {
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
	//calculate the size of the bounding box
	size = AABB_extent(ls);		
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
			m += LevelSetMethods::smoothed_heaviside(gr.ls, domain, i+domain.buffer, j+domain.buffer);
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
			c_x += LevelSetMethods::smoothed_heaviside(ls, domain, i+domain.buffer, j+domain.buffer)*(domain.X(i, j).x);
			c_y += LevelSetMethods::smoothed_heaviside(ls, domain, i+domain.buffer, j+domain.buffer)*(domain.X(i, j).y);
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
			inertia += LevelSetMethods::smoothed_heaviside(ls, domain, i+domain.buffer, j+domain.buffer)*(pow((i*domain.dx - c(0)),2) + pow((j*domain.dy - c(1)),2));
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
	LevelSet ls;
	ls.phi = matrix::Zero(domain.Nx+2*domain.buffer, domain.Ny+2*domain.buffer);

	//initialise the levelset with a large positive number
	for (int i=0; i<domain.Nx+2*domain.buffer; i++){
		for (int j=0; j<domain.Ny+2*domain.buffer; j++){
			ls.phi(i, j) = 1e6;
		}
	}

	for (int i=0; i<domain.Nx+2*domain.buffer; i++){
		for (int j=0; j<domain.Ny+2*domain.buffer; j++){
			double min_phi = 1e6;
			for (int a=0; a<static_cast<int>(particles.size()); a++){
				if (particles[a].ls.phi(i, j) < min_phi) min_phi = particles[a].ls.phi(i, j);
			}
			ls.phi(i, j) = min_phi;
		}
	}
	//LevelSetMethods::boundary_conditions(ls, domain);

	return ls;
}
////////////////////////////////////////////////
// THIS IS WRONG. 
//problem: vc*total time is not the total translation?
//Should store a displacement variable s and theta
//to track total displacement and rotation from initial axis
//the translation equation is integral of v(t) not final v.
//IDIOT :(
////////////////////////////////////////////////
LevelSet Particle::motion(const Domain2D& domain){
	//Temporary storage for levelset values at time t
	//The levelset can be updated using its initial values and current position
	LevelSet ls_t; //ls at time t
	ls_t.phi = matrix::Zero(domain.Nx+2*domain.buffer, domain.Ny+2*domain.buffer);

	//simultaneous translation and rotation
	//double pi = atan(1.0)*4;
	//double frequency = w/(2*pi);

	for(int i=0; i<domain.Nx+2*domain.buffer; i++){
		for (int j=0; j<domain.Ny+2*domain.buffer; j++){
			Coordinates originalpos(i*domain.dx - domain.buffer*domain.dx, j*domain.dy - domain.buffer*domain.dy);
			originalpos = LevelSetMethods::translation_reverse(originalpos, s); //translate the point
			originalpos = LevelSetMethods::rotation_reverse(originalpos, centroid, theta); //rotate the point
			ls_t.phi(i, j) = LevelSetMethods::interpolation_value(ls, domain, originalpos);
		}
		//std::cout << std::endl;
	}

	//move the nodes from ref, store a seperate node list
	for (int a=0; a<static_cast<int>(nodes.size()); a++){
		//std::cout << nodes[a].transpose() << std::endl;
		Coordinates node_a(ref_nodes[a]); //always start from reference nodes
		node_a = LevelSetMethods::rotation(node_a, centroid, theta);
		node_a = LevelSetMethods::translation(node_a, s); 
		nodes[a] = vector2(node_a.x, node_a.y);
	}
////////
	LevelSetMethods::fast_sweep(ls_t, domain);
	//LevelSetMethods::boundary_conditions(ls_t, domain);
	return ls_t;
}

vector2 Particle::cross(double w, const vector2& pos){
 	return vector2(-w*pos(1), w*pos(0));
}

double Particle::cross(const vector2& m, const vector2& n){
	return m(0)*n(1) - m(1)*n(0);
}

static double AABB_extent(const LevelSet& ls){
	double size=0;
	for (int i=0; i<ls.phi.rows(); i++){
		for (int j=0; j<ls.phi.cols(); j++){
			if (ls.phi(i, j) < size) size = ls.phi(i, j);
		}
	}
	//currently, size measure of the furthest point of the particle from its centroid.
	//double size to get the bounding box.
	size = -2*size;
	return size; 
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