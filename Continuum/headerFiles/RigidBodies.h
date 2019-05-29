#ifndef RIGIDBODIES_H_
#define RIGIDBODIES_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "Solvers.h"
#include "LevelSet.h"

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
	static void wall_collision(const vector2&, Particle&);
	static void newton_euler(Particle&, const Domain2D&, const vector2& torque, double force, double);
	static void initial_conditions(demTests&);
	static void solver(Moving_RB&, Domain2D&, double);
	static void output(const Moving_RB&, const Domain2D&);
	static void rigid_body_solver(demTests&, double);
	

};

#endif /* RIGIDBODIES_H_ */
