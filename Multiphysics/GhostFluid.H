#ifndef GHOSTFLUID_H_
#define GHOSTFLUID_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "RPsolvers/Solvers.h"
#include "RPsolvers/LevelSet.h"

//Note 1: Update is split from flux computation to allow identification of max wavespeed between both materials

//velocity to update levelset
/*
		if (sgn(phi(i)) < 0){
			velocity = U1(i, 1)/U1(i, 0); //U1 is the 'left' material
		}

		else{
			velocity = U2(i, 1)/U2(i, 0); //U2 (right) material is the real material
		}
*/

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
	void ghost_boundary_RP(); //RP based

	//HLLC
	//void initialise_HLLC(gfmTests);
	void initial_conditions_HLLC(EOS*, EOS*, gfmTests);
	void initial_conditions_HLLC(EOS*, EOS*, EOS*, gfmTests); //3 initial conditions, 2 discontinuities
	void update_levelset_HLLC(double);
	void solver_HLLC(EOS*, EOS*, gfmTests);
	void solver_HLLC(EOS*, EOS*, EOS*, gfmTests);
	void output_HLLC(EOS*, EOS*);
	void output_HLLC(EOS*, EOS*, EOS*);

	//MUSCL
	void initial_conditions_MUSCL(EOS*, EOS*, gfmTests);
	void initial_conditions_MUSCL(EOS*, EOS*, EOS*, gfmTests);
	void update_levelset_MUSCL(double);	
	void solver_MUSCL(EOS*, EOS*, gfmTests);
	void solver_MUSCL(EOS*, EOS*, EOS*, gfmTests);	
	void output_MUSCL(EOS*, EOS*);
	void output_MUSCL(EOS*, EOS*, EOS*);
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

};

/*
class GhostFluidMethods : public LevelSetFunction
{
	RPsolvers* var1;
	RPsolvers* var2;

	//domain wide dt and Smax after considering all materials
	double Smax;
	double dt;
	int count;

	public:
	//GhostFluidMethods(double, eulerTests, HLLC);
	//GhostFluidMethods(double, eulerTests, MUSCL);
	GhostFluidMethods(gfmTests);
	~GhostFluidMethods();
};
*/





#endif /* GHOSTFLUID_H_ */




