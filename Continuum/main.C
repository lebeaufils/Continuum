#include "headerFiles/Solvers.h"
#include "headerFiles/LevelSet.h"
#include "headerFiles/RigidBodies.h"
#include <typeinfo>


int main(void){


	demTests Tests(401);
	Tests.test9();

	RigidBodies::rigid_body_solver(Tests, 0.5);

 //error calc
	/*
	std::ofstream outfile;
	outfile.open("normal_error_square.txt");

	for (int i=1; i!=5; i++){
		int num = 40*pow(2, i);
		demTests Tests(num+1);
		Tests.test6();

		std::cout << "a" << std::endl;
		double h = 4.0/(num);
		vector2 n = Particle::normal_sum(Tests.var.particles.back().ls, Tests.domain);

		double magnitude = sqrt(n(1)*n(1) + n(0)*n(0));

		
		std::cout << h << '\t' << magnitude << std::endl;
		outfile << h << '\t' << num << '\t' << n.transpose() << '\t' << magnitude << std::endl;
	}
	//plot error in normal sum
	outfile.close();
	*/
}

