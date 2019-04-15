#ifndef RIGIDBODIES_H_
#define RIGIDBODIES_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "Solvers.h"

struct RigidBodies
{
	static void rigid_boundary(Euler1D, Euler1D&, int);
	static void initial_conditions();
	static void solver();
	static void output();

	static void rigid_body_solver();
};

#endif /* RIGIDBODIES_H_ */
