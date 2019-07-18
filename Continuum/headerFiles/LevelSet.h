#ifndef LEVELSET_H_
#define LEVELSET_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <map>

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
	static Coordinates translation(const Coordinates&, const vector2&);
	static Coordinates rotation(const Coordinates&, const vector2&, double);
	static Coordinates translation_reverse(const Coordinates&, const vector2&);
	static Coordinates rotation_reverse(const Coordinates&, const vector2&, double);
	static LevelSet motion(const LevelSet&, const Domain2D&, const vector2&, const vector2&, double);

};

#endif /* LEVELSET_H_ */