#ifndef RIGIDBODIES_H_
#define RIGIDBODIES_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "Solvers.h"

struct RigidBodies
{
	static void fast_sweep(const LevelSet&, RB_2D&, const Domain2D&, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic>&);
	static vector4 reflected_state(const RB_2D&, int, int, const vector2&);
	static void boundary_conditions(RB_2D&, const Domain2D&);
	static void initial_conditions(rigidTests&);
	static void solver();
	static void output();

	static void rigid_body_solver();
};

#endif /* RIGIDBODIES_H_ */
