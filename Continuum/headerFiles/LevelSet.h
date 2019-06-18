#ifndef LEVELSET_H_
#define LEVELSET_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

//#include "eulerTests.h"
#include "Variables.h"

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;


class LevelSetMethods
{
public:

	static int get_sgn(double);
	
	//-----------------------------------------------------
	//1-Dimensional
	//-----------------------------------------------------
	static void initialise(LevelSet&, const Domain1D&);

	static void boundary_conditions(LevelSet&, const Domain1D&);
	static void signed_distance_function(LevelSet&, const Domain1D&, double);
	static void signed_distance_function(LevelSet&, const Domain1D&, double, double, int); //2 discontinuties
	//void signed_distance_function_1D_3(int); //3 discontinuities
	
	//first order upwind hamilton-jacobi method
	//static double HJ_FirstOrder(double, double, int); //velocity
	//void HJ_ENO(VectorN2, double);
	//void HJ_WENO(VectorN2, double);

	//reinitialisation
	//static void reinitialisation();
	//void reconstruction();

	//-----------------------------------------------------
	//2-Dimensional
	//-----------------------------------------------------
	static void boundary_conditions(LevelSet&, const Domain2D&);
	static void initialise(LevelSet&, const Domain2D&, const Polygon&);
	static void initialise_circle(LevelSet&, const Domain2D&, double, double, double);

	static void fast_sweep(LevelSet&, const Domain2D&);

	static vector2 normal(const LevelSet&, const Domain2D&, int, int);
	static double interpolation_value(const LevelSet&, const Domain2D&, const Coordinates&);
	static vector2 interpolation_gradient(const LevelSet&, const Domain2D&, const Coordinates&);

	//-----------------------------------------------------
	//Utility
	//-----------------------------------------------------
	static double smoothed_heaviside(const LevelSet&, const Domain2D&, int, int);
	static double smoothed_delta(const LevelSet&, const Domain2D&, int, int);
	//binary operators
	static LevelSet merge(const std::vector<LevelSet>&, const Domain2D&); //union of levelsets
	//static LevelSet intersection(const std::vector<LevelSet>&);
	//static LevelSet difference(const std::vector<LevelSet>&);

	//-----------------------------------------------------
	//Forces
	//-----------------------------------------------------
	static vector2 force(const Euler2D&, const LevelSet&, const Domain2D&);
	static double torque (const Euler2D&, const LevelSet&, const Domain2D&, const vector2&);

	//-----------------------------------------------------
	//Motion
	//-----------------------------------------------------
	static Coordinates translation(const Coordinates&, const vector2&, double);
	static Coordinates rotation(const Coordinates&, const vector2&, double, double);
	static Coordinates translation_reverse(const Coordinates&, const vector2&, double);
	static Coordinates rotation_reverse(const Coordinates&, const vector2&, double, double);
	static LevelSet motion(const LevelSet&, const Domain2D&, const vector2&, const vector2&, double, double);

};

//------------------------------------------------------------
//Particles 
//------------------------------------------------------------
//Particles are represented by a levelset function
struct Particle //NEEDS WORK
{
	//RB_2D should store particles rather than levelsets.
	//each particle can have its reference levelset and nodes
	double density = 3000;

	LevelSet ls;
	vector2 centroid; //initial, fixed
	vector2 centre;

	vector2 vc; //translational velocity
	double w; //angular velocity
	std::vector<vector2> nodes;
	std::vector<vector2> ref_nodes;
	//???If the mass of the Particle is nott important

	//Stiffness
	double miu; //interparticle friction coefficient
	double k_n; //normal contact stiffness
	double k_s; //shear contact stiffness
	//Forces
	vector2 force; //Accumulator for linear force
	double torque; //Accumuluator for torque

	Particle() : ls(), centroid(0, 0), centre(0, 0), vc(0, 0), w(0), nodes(0), ref_nodes(0), miu(0.26), k_n(1000), k_s(1000), force(0, 0), torque(0) {}
	Particle(const Domain2D&, const Coordinates&, double);
	Particle(const Polygon&, const Domain2D&);
	//Particle(const Particle&) //copy constructor
	Particle(const Particle& gr);
	~Particle() {};

	//void initialise(const Polygon&, const Domain2D&);
	//Particle(double d); //no real need to specify density for rigid bodies
	double diameter(); //estimated particle diameter, equidiv?
	//void check_contact(); //?
	void set_velocity(const vector2&, double);
	LevelSet motion(const Domain2D&, double);

	//-----------------------------------------------------
	//Inertial properties
	//-----------------------------------------------------
	static double mass(const Particle&, const Domain2D&);
	static vector2 center_of_mass(const LevelSet&, const Particle&, const Domain2D&);
	static double moment_of_inertia(const LevelSet&, const Particle&, const Domain2D&);
	static vector2 velocity(const Coordinates&, const Particle&); //vb

	static vector2 cross(double, const vector2&);
	static double cross(const vector2&, const vector2&);
	static LevelSet merge(const std::vector<Particle>&, const Domain2D&);
};

struct Moving_RB
{
	//Rigid body system of moving particles
	//collection of "particles" which contain their own levelset and discretised nodes
	std::vector<Particle> particles;
	Euler2D fluid;

	LevelSet combinedls;
	//Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal;

	Moving_RB() : particles(0), fluid(), combinedls() {} //, normal(0, 0) {}
	~Moving_RB() {};

	void add_sphere(const Domain2D&, const Coordinates&, double); //center and radius
	void add_particle(const Polygon&, const Domain2D&);

	int getNBodies(); //returns number of rigid body particles

private:
	Moving_RB(const Moving_RB& rbsystem) : particles(rbsystem.particles), fluid(rbsystem.fluid), combinedls(rbsystem.combinedls) {}

};

#endif /* LEVELSET_H_ */