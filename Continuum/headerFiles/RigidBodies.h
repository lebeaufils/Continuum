#ifndef RIGIDBODIES_H_
#define RIGIDBODIES_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "Solvers.h"
#include "LevelSet.h"
#include <string>
#include <boost/math/constants/constants.hpp>

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
	static vector2 compute_tangential_force(vector2&, Particle&, double, double, const vector2&, const vector2&);
	static double compute_normal_force(const Particle&, double, const vector2&, double);
	static vector2 compute_tangential_force(vector2&, Particle&, double, double, const vector2&, const vector2&, double);
	static double compute_torque(const Particle&, const vector2&, const vector2&, const vector2&);

	//static void contact_detection(const Domain2D&, std::vector<Particle>&, const std::vector<LevelSet>&, double); //node to levelset contact check
	//need to update spring implementation
	static std::unordered_map<int, vector2>::iterator find_spring(std::unordered_map<int, vector2>, int);
	static void delete_spring(std::unordered_map<int, vector2>, int);
	static void contact_detection(const Domain2D&, Moving_RB&, double);
	static void wall_collision(vector2&, Particle&, const vector2&, double, const vector2&, double);
	static void particle_collision(vector2&, Particle&, Particle&, const vector2&, double, const vector2&, double);
	//
	static void fluid_forces(const Domain2D&, const Euler2D&, std::vector<Particle>&);
	static void fluid_forces(const Domain2D&, const Euler2D&, std::vector<Particle>&, double); //with gravity
	static void newton_euler(Particle&, const Domain2D&, const vector2& torque, double force, double);
	static void update_displacements(Particle&, const Domain2D&, Moving_RB&, double);
	static void update_displacements_forced_motion(Particle&, const Domain2D&, Moving_RB&, double);
	static void initial_conditions(demTests&);
	static void subcycling(Moving_RB&, const Domain2D&, const Euler2D&, double, double);
	static void forced_motion(Moving_RB&, const Domain2D&, const Euler2D&, double);
	static void solver(Moving_RB&, Domain2D&, double);
	static void output(const Moving_RB&, const Domain2D&, std::string, std::string);
	static void output_levelset(const Moving_RB&, const Domain2D&, std::string, std::string);
	static void rigid_body_solver(demTests&, double);
};

#endif /* RIGIDBODIES_H_ */
