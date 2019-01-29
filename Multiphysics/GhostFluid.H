#ifndef GHOSTFLUID_H_
#define GHOSTFLUID_H_

#include "Solvers.h"
#include "LevelSet.h"

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

/*class GhostFluidMethods : public HLLC, public MUSCL, public LevelSetFunction
{

	//Accessing members of HLLC/MUSCL as such
	//A::word = a;
	//B::word = b;
public:
	GFM(double, eulerTests);

	void initial_conditions(EOS*, EOS*, eulerTests); //EOS is polymorphic
	void ghost_boundary(); //original GFM
	void ghost_boundary_RP(); //RP based
	void update_levelset();
	//void update_levelset_ENO();
	//void update_levelset_WENO();
	//void update_levelset_TVD();
	void solver();
};
*/
class GhostFluidMethods : public LevelSetFunction
{
	RPsolvers* var1;
	RPsolvers* var2;

	//domain wide dt and Smax after considering all materials
	double Smax;
	double dt;
public:
	//GhostFluidMethods(double, eulerTests, HLLC);
	//GhostFluidMethods(double, eulerTests, MUSCL);
	GhostFluidMethods(double, gfmTests);
	~GhostFluidMethods();

	void ghost_boundary(gfmTests, matrix, matrix, int); //original GFM
	void ghost_boundary_RP(); //RP based

	template <typename T> void initial_conditions(T, EOS*, EOS*, gfmTests){static_assert(sizeof(T) == 0, "Specialised calls only")}
	template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, EOS*, gfmTests);
	void update_levelset(double);
	//void update_levelset_ENO();
	//void update_levelset_WENO();
	//void update_levelset_TVD();
	void solver();
}

#endif /* GHOSTFLUID_H_ */