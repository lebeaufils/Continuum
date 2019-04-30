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
	static void initialise(LevelSet&, const Domain2D&, Polygon&);
	static void initialise_circle(LevelSet&, Domain2D, double, double, double);

	static void fast_sweep(LevelSet&, Domain2D);
	static vector2 normal(LevelSet, Domain2D, int, int);
};


#endif /* LEVELSET_H_ */