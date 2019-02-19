#ifndef GHOSTFLUID_H_
#define GHOSTFLUID_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "RPsolvers/Solvers.h"
#include "RPsolvers/LevelSet.h"


class GhostFluidMethods : public LevelSetFunction
{
	MUSCL* var1;
	MUSCL* var2;

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
		delete var1;
		delete var2;
	};

	void ghost_test_rightghost(MUSCL*, EOS*, MUSCL*, EOS*, int);
	void ghost_test_leftghost(MUSCL*, EOS*, MUSCL*, EOS*, int);
	
	void ghost_boundary(MUSCL*, EOS*, MUSCL*, EOS*, int); //original GFM
	void ghost_boundary_RP(MUSCL*, EOS*, MUSCL*, EOS*, int); //RP based

	//MUSCL
	void update_levelset(double);	
	//Initial conditions
	void initial_conditions(EOS*, EOS*, gfmTests);
	void initial_conditions(EOS*, EOS*, EOS*, gfmTests);
	void initial_conditions(EOS*, EOS*, EOS*, EOS*, gfmTests);
	//Solver
	void solver(gfmTests);
		//void solver(EOS*, EOS*, gfmTests);
		//void solver(EOS*, EOS*, EOS*, gfmTests);	
		//void solver(EOS*, EOS*, EOS*, EOS*, gfmTests);

	//output
	void output(EOS*, EOS*);
	void output(EOS*, EOS*, EOS*);
	void output(EOS*, EOS*, EOS*, EOS*);

	//RP based GFM solver
	//initial conditions
	void initial_conditions_RP(EOS*, EOS*, gfmTests);
	void initial_conditions_RP(EOS*, EOS*, EOS*, gfmTests);
	void initial_conditions_RP(EOS*, EOS*, EOS*, EOS*, gfmTests);
	//solver
	void solver_RP(gfmTests);	
	//void solver_RP(EOS*, EOS*, gfmTests);
	//void solver_RP(EOS*, EOS*, EOS*, EOS*, gfmTests);

	//void update_levelset_ENO();
	//void update_levelset_WENO();
	//void update_levelset_TVD();

	//void reinitialise();

	//EXACT
	double f(double, EOS*, EOS*);
	double fprime(double, EOS*, EOS*);
	void check_pressure_pos_condition(EOS*, EOS*);
	double newton_raphson(double, EOS*, EOS*);
	double compute_star_pressure(EOS*, EOS*);
	double compute_star_velocity(double, EOS*, EOS*);
	double compute_shock_density(double, EOS*);
	double compute_rarefraction_density(double, EOS*);
	//void exact_solver(gfmTests, EOS*, EOS*); 
	//void exact_solver(gfmTests, EOS*, EOS*, EOS*);
	void exact_solver(gfmTests);

	double compute_star_pressure_SG(StiffenedGas, StiffenedGas);
	void exact_solver_SG(gfmTests); 
	double compute_shock_density_SG(double, StiffenedGas);
	double compute_rarefraction_density_SG(double, StiffenedGas);

	void ghost_boundary_RP_SG(MUSCL*, StiffenedGas, MUSCL*, StiffenedGas, int); //RP based
	void initial_conditions_RP_SG(StiffenedGas&, StiffenedGas&, gfmTests);
	void solver_RP_SG(gfmTests);

	//void ghost_boundary_RP_SG(MUSCL*, StiffenedGas*, MUSCL*, StiffenedGas*, int); //RP based

};






#endif /* GHOSTFLUID_H_ */




