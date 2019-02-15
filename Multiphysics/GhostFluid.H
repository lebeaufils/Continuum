#ifndef GHOSTFLUID_H_
#define GHOSTFLUID_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "RPsolvers/Solvers.h"
#include "RPsolvers/LevelSet.h"


class GhostFluidMethods : public LevelSetFunction
{
	RPsolvers* var1;
	RPsolvers* var2;

	//domain wide dt and Smax after considering all materials
	double CFL;
	double Smax;
	double dt;
	int count;

	//constants storage
	matrix C;
public:
	//GhostFluidMethods(double, eulerTests, HLLC);
	//GhostFluidMethods(double, eulerTests, MUSCL);
	GhostFluidMethods(double c, gfmTests);
	~GhostFluidMethods() {
		delete var1; delete var2;
	};

	
	void ghost_boundary(RPsolvers*, EOS*, RPsolvers*, EOS*, int); //original GFM
	void ghost_boundary_RP(RPsolvers*, EOS*, int, RPsolvers*, EOS*, int, int); //RP based

	//HLLC
	//void initialise_HLLC(gfmTests);
	void initial_conditions_HLLC(EOS*, EOS*, gfmTests);
	void initial_conditions_HLLC(EOS*, EOS*, EOS*, gfmTests); //3 initial conditions, 2 discontinuities
	//void initial_conditions_HLLC(EOS*, EOS*, EOS*, EOS*, gfmTests); 
	void update_levelset_HLLC(double);
	void solver_HLLC(EOS*, EOS*, gfmTests);
	void solver_HLLC(EOS*, EOS*, EOS*, gfmTests);
	//void solver_HLLC(EOS*, EOS*, EOS*, EOS*, gfmTests); //too diffusive
	void output_HLLC(EOS*, EOS*);
	void output_HLLC(EOS*, EOS*, EOS*);
	//void output_HLLC(EOS*, EOS*, EOS*, EOS*);

	//MUSCL
	void initial_conditions_MUSCL(EOS*, EOS*, gfmTests);
	void initial_conditions_MUSCL(EOS*, EOS*, EOS*, gfmTests);
	void initial_conditions_MUSCL(EOS*, EOS*, EOS*, EOS*, gfmTests);
	void update_levelset_MUSCL(double);	
	void solver_MUSCL(EOS*, EOS*, gfmTests);
	void solver_MUSCL(EOS*, EOS*, EOS*, gfmTests);	
	void solver_MUSCL(EOS*, EOS*, EOS*, EOS*, gfmTests);

	//RP based GFM solver
	void initial_conditions_MUSCL_RP(EOS*, EOS*, gfmTests);
	void solver_MUSCL_RP(EOS*, EOS*, gfmTests);

	//output
	void output_MUSCL(EOS*, EOS*);
	void output_MUSCL(EOS*, EOS*, EOS*);
	void output_MUSCL(EOS*, EOS*, EOS*, EOS*);

	//EXACT
	void y_constants(gfmTests);
	double fk(double, vector, int);
	double fprimek(double, vector, int);
	double f(double, vector, int, vector, int);
	double fprime(double, vector, int, vector, int);
	void check_pressure_pos_condition(vector, int, vector, int);
	double newton_raphson(double, vector, int, vector, int);
	double compute_star_pressure(vector, int, vector, int);
	double compute_star_velocity(double, vector, int, vector, int);
	double compute_shock_density(double, vector, int);
	double compute_rarefraction_density(double, vector, int);
	void exact_solver(gfmTests); 
	//remember to check that internal energy is divided by density (specific)

	//template <typename T> void initial_conditions(T, EOS*, EOS*, gfmTests){static_assert(sizeof(T) == 0, "Specialised calls only")}


	//template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, EOS*, gfmTests);

	//void update_levelset_ENO();
	//void update_levelset_WENO();
	//void update_levelset_TVD();

	//void reinitialise();


	//TESTING
	double testingf(double, EOS*, EOS*);
	double testingfprime(double, EOS*, EOS*);
	void testingcheck_pressure_pos_condition(EOS*, EOS*);
	double testingnewton_raphson(double, EOS*, EOS*);
	double testingcompute_star_pressure(EOS*, EOS*);
	double testingcompute_star_velocity(double, EOS*, EOS*);
	double testingcompute_shock_density(double, EOS*);
	double testingcompute_rarefraction_density(double, EOS*);
	void testingexact_solver(gfmTests, EOS*, EOS*); 

};






#endif /* GHOSTFLUID_H_ */




