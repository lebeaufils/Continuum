#ifndef GHOSTFLUID_H_
#define GHOSTFLUID_H_

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
public:
	//GhostFluidMethods(double, eulerTests, HLLC);
	//GhostFluidMethods(double, eulerTests, MUSCL);
	GhostFluidMethods(double c, gfmTests);
	~GhostFluidMethods() {
		delete var1; delete var2;
	};

	void initialise_HLLC(gfmTests);
	//void initialise_EXACT():
	//void initialise_MUSCL();

	void ghost_boundary(RPsolvers*, EOS*, RPsolvers*, EOS*, int); //original GFM
	void ghost_boundary_RP(); //RP based

	//template <typename T> void initial_conditions(T, EOS*, EOS*, gfmTests){static_assert(sizeof(T) == 0, "Specialised calls only")}
	void initial_conditions_HLLC(EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<HLLC>(HLLC, HLLC, EOS*, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, EOS*, gfmTests);
	//template <> void initial_conditions<MUSCL>(MUSCL, MUSCL, EOS*, EOS*, EOS*, gfmTests);
	void update_levelset(double);
	//void update_levelset_ENO();
	//void update_levelset_WENO();
	//void update_levelset_TVD();
	void solver(EOS*, EOS*, gfmTests);
	void output(EOS*, EOS*);


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




