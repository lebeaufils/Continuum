#include "headerFiles/Solvers.h"
#include "headerFiles/LevelSet.h"
#include "headerFiles/RigidBodies.h"
#include <typeinfo>


int main(void){


	demTests Tests(401);
	Tests.test7();

	RigidBodies::rigid_body_solver(Tests, 0.5);

/* //error calc
	std::ofstream outfile;
	outfile.open("normal_error.txt");

	for (int i=1; i!=20; i++){
		demTests Tests(40*i+1);
		Tests.test6();

		double h = 4.0/(40*i);
		vector2 n = Particle::normal_sum(Tests.var.particles.back().ls, Tests.domain);

		double magnitude = sqrt(n(1)*n(1) + n(0)*n(0));

		
		std::cout << h << '\t' << magnitude << std::endl;
		outfile << 40*i << '\t' << n.transpose() << '\t' << magnitude << std::endl;
	}
	//plot error in normal sum
	outfile.close();
*/
}

