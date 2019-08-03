#include "headerFiles/Solvers.h"
#include "headerFiles/LevelSet.h"
#include "headerFiles/RigidBodies.h"
#include <typeinfo>


int main(void){


	demTests Tests(201);
	Tests.test4();

	//RigidBodies::initial_conditions(Tests);
	RigidBodies::rigid_body_solver(Tests, 0.5);

/*	Polygon poly;
	Coordinates center(0.5,0.5);
	try{
		poly.create_square(Tests.domain, 0.3, center);
	}
	catch (const char c){
		std::cout << c << std::endl;
	}
	Particle gr(poly, Tests.domain);

	RigidBodies::newton_euler(Tests.var.fluid, gr, Tests.domain, 0.1);
	std::cout << gr.vc.transpose() << std::endl;
	std::cout << gr.w << std::endl;
*/

	/*vector2 c =  LevelSetMethods::center_of_mass(Tests.var.levelsets[0], Tests.domain, poly);
	std::cout << "Center of mass = " << c.transpose() << std::endl;
	for (int i=0; i<static_cast<int>(poly.vertices.size()); i++){
		vector2 v(poly.vertices[i].x, poly.vertices[i].y);
		vector2 r = Rotor2::rotate_about(v, c, 1, 4*atan(1.)/3);
		poly.vertices[i].x = r(0);
		poly.vertices[i].y = r(1);
	}

	poly.output(Tests.domain);
	std::cout << poly.vertices[0].x << '\t' << poly.vertices[0].y << std::endl;
	*/
}

