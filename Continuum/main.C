#include "headerFiles/Solvers.h"
#include <typeinfo>


int main(void){

	//eulerTests2D Tests(101);
	//Tests.test3();
	//MUSCL::muscl_solver(Tests, 0.9);

	//eulerTests Test1d(101);
	//Test1d.test1_stationary();
	//Tests.switch_resolution();
	//Tests.switch_test();

	//EXACT theExact(Test1d);
	//theExact.solver(Test1d);

	rigidTests Tests(11);
	Tests.test1();

	Coordinates P1(0.5, 0.8);
	Coordinates P2(0.8, 0.25);

	std::vector<Pos_Index> indices;

	indices = Bresenham::line_algorithm(Tests.domain, P1, P2);

	std::cout << std::endl;
	std::cout << std::endl;

	//for (int i=0; i<static_cast<int>(indices.size()); i++){
		//Tests.domain.X(indices[i].i, indices[i].j).display();
		//std::cout << indices[i].i << ", " << indices[i].j << '\t';
	//}
	for (int i=0; i<static_cast<int>(indices.size()); i++){
		Tests.domain.X(indices[i].i, indices[i].j).display();
	}
}

