#ifndef RIGIDBODIES_H_
#define RIGIDBODIES_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "Solvers.h"
#include "LevelSet.h"
#include <string>

struct RigidBodies
{

	static void boundary_conditions(vecarray&, const Domain2D&);
	//---------------------------------------------------------------
	// Stationary rigidbody
	//---------------------------------------------------------------
	static void fast_sweep(const LevelSet&, Stationary_RB&, const Domain2D&, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic>&);
	static void reflected_state(Stationary_RB&, const vecarray&, int i, int j, const vector2&);
	static void initial_conditions(rigidTests&);
	static void solver(Stationary_RB&, Domain2D&, double);
	static void output(const Stationary_RB&, const Domain2D&);

	static void rigid_body_solver(rigidTests&, double);

	//---------------------------------------------------------------
	// Moving rigidbody and Forces
	//---------------------------------------------------------------
	//static void boundary_conditions(Stationary_RB&, const Domain2D&);
	//treat collision with domain boundaries as collisions?
	static void fast_sweep(const LevelSet&, const Particle&, Moving_RB&, const Domain2D&, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic>&);
	static void reflected_state(Moving_RB&, const vecarray&, int i, int j, const vector2&, const vector2&);
	//collisions
	static double compute_normal_force(const Particle&, double, const vector2&);
	static vector2 compute_tangential_force(int, Particle&, double, double, double, const vector2&, const vector2&);
	static double compute_torque(const Particle&, const vector2&, const vector2&, const vector2&);

	static void contact_detection(const Domain2D&, std::vector<Particle>&, const std::vector<LevelSet>&, double); //node to levelset contact check

	static void wall_collision(Particle&, const vector2&, double, const vector2&, double);
	static void particle_collision(int j, Particle&, Particle&, const vector2&, double, const vector2&, double);
	//
	static void newton_euler(const LevelSet&, Particle&, const Domain2D&, const vector2& torque, double force, double);
	static void update_displacements(Particle&, double);
	static void initial_conditions(demTests&);
	static void subcycling(Moving_RB&, const Domain2D&, std::vector<LevelSet>&, double, double);
	static void solver(Moving_RB&, Domain2D&, double);
	static void output(const Moving_RB&, const Domain2D&, std::string, std::string);
	static void output_levelset(const Moving_RB&, const Domain2D&, std::string, std::string);
	static void rigid_body_solver(demTests&, double);
};

#endif /* RIGIDBODIES_H_ */
