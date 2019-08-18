#include "../headerFiles/eulerTests.h"

int standardTests::Get_Switch(){
	int switch_value;
	do{
		std::cin >> switch_value;

		bool find;
		if (switch_value > 0 && switch_value < 10) find = 1;
		else find = 0;

		if (find == 1){
			return switch_value;
			break;
		}
		else{ //error if find fails
			std::cout << "Invalid input." << std::endl;
			std::cin.clear();
			std::cin.ignore();
			return Get_Switch();
		}

	}while(true);
}

void standardTests::switch_resolution(Domain1D &domain){ //1D
	std::cout << "Resolution options" << std::endl
		<< "1. Low (100 Cells)" << std::endl
		<< "2. Mediumn (200 Cells)" << std::endl
		<< "3. High (400 Cells)" << std::endl
		<< "4. Very High (1000 Cells)" << std::endl
		<< "5. Exit" << std::endl;

	int a = Get_Switch();

	switch(a){
		case 0:
		case 1:
			domain.N = 100;
			break;
		case 2:
			domain.N = 200;
			break;
		case 3:
			domain.N = 400;
			break;
		case 4:
			domain.N = 1000;
			break;
		case 5:
		case 6:
		case 7:
		case 8:
		case 9:
			exit(0);
	}
}

//-----------------------------------------

void eulerTests::test1(){ //sod's tube test
	vector Left(1.0, 0.75, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.3;
	domain.tstop = 0.2;
	domain.dx = domain.L/(domain.N-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test1_stationary(){
	vector Left(1.0, 0.0, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	domain.tstop = 0.25;
	domain.dx = domain.L/(domain.N-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test2(){ //Two rarefraction waves, vaccum test
	vector Left(1.0, -2.0, 0.4);
	vector Right(1.0, 2.0, 0.4);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	domain.tstop = 0.15;
	domain.dx = domain.L/(domain.N-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test3(){ //strong shock wave of shock Mach number 198
	vector Left(1.0, 0.0, 1000.0);
	vector Right(1.0, 0.0, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	domain.tstop = 0.012;
	domain.dx = domain.L/(domain.N-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test4(){ //Three strong discontinuities travelling to the right.
	vector Left(5.99924, 19.5975, 460.894);
	vector Right(5.99242, -6.19633, 46.0950);

	initialL = Left;
	initialR = Right;

	x0 = 0.4;
	domain.tstop = 0.035;
	domain.dx = domain.L/(domain.N-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test5(){ //slowly moving contact discontinuities
	vector Left(1.0, -19.59745, 1000.0);
	vector Right(1.0, -19.59745, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.8;
	domain.tstop = 0.012;
	domain.dx = domain.L/(domain.N-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::switch_test(){
	std::cout << "Test options" << std::endl
		<< "Toro's Shock Tube tests" << std::endl
		<< "1. Test 1" << std::endl
		<< "2. Test 2" << std::endl
		<< "3. Test 3" << std::endl
		<< "4. Test 4" << std::endl
		<< "5. Test 5" << std::endl
		<< "6. Exit" << std::endl;

	int a = Get_Switch();

	switch(a){
		case 0:
		case 1:
			test1();
			break;
		case 2:
			test2();
			break;
		case 3:
			test3();
			break;
		case 4:
			test4();
			break;
		case 5:
			test5();
			break;
		case 6:
		case 7:
		case 8:
		case 9:
			exit(0);
	}
}


//--------------------------------------------------------------------------------------------
//	2-D problems
//--------------------------------------------------------------------------------------------
void eulerTests2D::test1(){ //sod's tube test, x-aligned
	vector4 Left(1.0, 0.0, 0.0, 1.0); //density, velocity, pressure
	vector4 Right(0.125, 0.0, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	//x0 = 0.5;
	domain.tstop = 0.25;

	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;

	for (int i=0; i<domain.Nx; i++){
		x = i*domain.dx;
		for (int j=0; j<domain.Ny; j++){
			y = j*domain.dy;
			if (x <= 0.5){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
		}
	}
	//std::cout << interface << std::endl;
}

void eulerTests2D::test2(){ //sod's tube test, x-aligned
	vector4 Left(1.0, 0.0, 0.0, 1.0); //density, velocity, pressure
	vector4 Right(0.125, 0.0, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	//x0 = 0.5;
	domain.tstop = 0.25;

	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;

	for (int i=0; i<domain.Nx; i++){
		x = i*domain.dx;
		for (int j=0; j<domain.Ny; j++){
			y = j*domain.dy;
			if (y <= 0.5){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
		}
	}
	//std::cout << interface << std::endl;
}

void eulerTests2D::test3(){ //sod's tube test, x-aligned
	vector4 Left(1.0, 0.0, 0.0, 1.0); //density, velocity, pressure
	vector4 Right(0.125, 0.0, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	//x0 = 0.5;
	domain.tstop = 0.25;

	domain.Lx = domain.Ly = sqrt(0.5);

	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;
	for (int i=0; i<domain.Nx; i++){
		x = i*domain.dx;
		for (int j=0; j<domain.Ny; j++){
			y = j*domain.dy;
			if ((x + y) <= sqrt(0.5)){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
			//note the (i, j) confusion here
			//x sweep = constant row, changing columns
		}
	}
	//std::cout << interface << std::endl;
}

void eulerTests2D::test4(){ //sod's tube test, x-aligned
	vector4 in(1.0, 0.0, 0.0, 1.0); //density, velocity, pressure
	vector4 out(0.125, 0.0, 0.0, 0.1);

	initialL = out;
	initialR = in;

	//x0 = 0.5;
	domain.tstop = 0.25;

	domain.Lx = domain.Ly = 2.0;

	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);

	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;
	for (int i=0; i<domain.Nx; i++){
		x = i*domain.dx;
		for (int j=0; j<domain.Ny; j++){
			y = j*domain.dy;
			if ((pow((x-1),2) + pow((y-1),2)) >= pow(0.4, 2)){
				interface(i, j) = false; //outside
			}
			else {
				interface(i, j) = true;
			}
			//note the (i, j) confusion here
			//x sweep = constant row, changing columns
		}
	}
	//std::cout << interface << std::endl;
}

//--------------------------------------------------------------------------------------------
//	Rigid Bpdy Tests
//--------------------------------------------------------------------------------------------
void rigidTests::test1(){

	//Shock location
	//double x_s = 0.3;

	//------------------------------------------------------------
	//initial conditions for shock and unshocked fluid
	vector4 shocked(1.3764, 0.394, 0.0, 1.5698); //density, velocity, pressure
	vector4 unshocked(1.0, 0.0, 0.0, 1.0);

	initialL = shocked;
	initialR = unshocked;

	domain.tstop = 0.25;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (domain.X(i, j).x <= 0.2){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	var.add_levelset();
	//generating the level set as a circle
	LevelSetMethods::initialise_circle(var.levelsets[0], domain, 0.5, 0.5, 0.2);

	//domain.display_grid();

	//std::cout << std::endl;

	//var.levelsets[0].display_grid();

	//a single interface between rigid body and fluid
	//essentially, this reduces the computational domain by bringing forward the boundary
	//std::cout << interface << std::endl;
}

void rigidTests::test2(){

	//Shock location
	//double x_s = 0.3;

	//------------------------------------------------------------
	//initial conditions for shock and unshocked fluid
	vector4 shocked(1.3764, 0.394, 0.0, 1.5698); //density, velocity, pressure
	vector4 unshocked(1.0, 0.0, 0.0, 1.0);

	initialL = shocked;
	initialR = unshocked;

	domain.tstop = 0.25;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (domain.X(i, j).x <= 0.2){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	Polygon poly;
	Coordinates center(0.5,0.5);
	try{
		poly.create_square(domain, 0.3, center);
	}
	catch (const char c){
		std::cout << c << std::endl;
	}

	var.add_levelset();
	//generating the level set as a circle
	LevelSetMethods::initialise(var.levelsets[0], domain, poly);

}

void rigidTests::test3(){

	//Shock location
	//double x_s = 0.3;

	//------------------------------------------------------------
	//initial conditions for shock and unshocked fluid
	vector4 shocked(1.3764, 0.394, 0.0, 1.5698); //density, velocity, pressure
	vector4 unshocked(1.0, 0.0, 0.0, 1.0);

	initialL = shocked;
	initialR = unshocked;

	domain.tstop = 0.25;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (domain.X(i, j).x <= 0.2){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	Polygon poly;
	try{
		poly.create(domain, 0.3, 10);
	}
	catch (const char c){
		std::cout << c << std::endl;
	}

	var.add_levelset();
	//generating the level set as a circle
	LevelSetMethods::initialise(var.levelsets[0], domain, poly);

}

void rigidTests::test4(){

	//Shock location
	//double x_s = 0.3;

	//------------------------------------------------------------
	//initial conditions for shock and unshocked fluid
	vector4 shocked(1.3764, 0.394, 0.0, 1.5698); //density, velocity, pressure
	vector4 unshocked(1.0, 0.0, 0.0, 1.0);

	initialL = shocked;
	initialR = unshocked;

	domain.tstop = 0.25;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (domain.X(i, j).x <= 0.2){
				interface(i, j) = false;
			}
			else {
				interface(i, j) = true;
			}
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	Polygon poly(1.0);
	try{
		poly.create_from_file(domain);
	}
	catch (const char c){
		std::cout << c << std::endl;
	}

	var.add_levelset();
	//generating the level set as a circle
	LevelSetMethods::initialise(var.levelsets[0], domain, poly);
	//std::cout << "mass = " << LevelSetMethods::mass(var.levelsets[0], domain, poly) << std::endl;
	//std::cout << "Xc = " << LevelSetMethods::center_of_mass(var.levelsets[0], domain, poly).transpose() << std::endl;
	//std::cout << "moment of inertia = " << LevelSetMethods::moment_of_inertia(var.levelsets[0], domain, poly) << std::endl;
}

//--------------------------------------------------------------------------------------------
//	Moving Particles (DEM) Tests
//--------------------------------------------------------------------------------------------
void demTests::test1(){

	//Vaccuum

	vector4 unshocked(0.0, 0.0, 0.0, 0.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 0.2;
	domain.Lx = 4.0; domain.Ly = 4.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	//var.add_sphere(domain, Coordinates(0.2,0.5), 0.1);
	//var.particles.back().set_velocity(vector2(2, 0), 0);
	//var.add_sphere(domain, Coordinates(0.8,0.5), 0.1);
	//var.particles.back().set_velocity(vector2(-2, 0), 0);
	Polygon poly;
	poly.create_from_file(domain,"vertices_zero.txt", vector2(1.5, 2.0));
	var.add_particle(poly, domain);
	var.particles.back().set_velocity(vector2(4, 0), 0);

	Polygon poly1;
	poly1.create_from_file(domain,"vertices_zero.txt", vector2(2.5, 2.0));
	var.add_particle(poly1, domain);
	var.particles.back().set_velocity(vector2(-4, 0), 0);


	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 100;//1e6;
		i->k_s = 100;//1e6;
		i->k_c = 100;//1e6;
		i->density = 0.625; //3000;
	}
}

void demTests::test2(){

	//Vaccuum
	vector4 unshocked(0.0, 0.0, 0.0, 0.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 0.5;
	domain.Lx = 4.0; domain.Ly = 4.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------

	var.add_sphere(domain, Coordinates(0.5,0.5), 0.5);
	var.add_sphere(domain, Coordinates(1.5, 0.5), 0.5);
	var.add_sphere(domain, Coordinates(2.5, 0.5), 0.5);
	var.add_sphere(domain, Coordinates(3.5, 0.5), 0.5);

	var.gravity = 10.0;
	var.add_sphere(domain, Coordinates(2.0, 2.0), 0.5);
	var.particles.back().set_velocity(vector2(0.0, -2.0), 0.0);
	//var.add_sphere(domain, Coordinates(1.2, 0.3), 0.2);
	//var.add_sphere(domain, Coordinates(1.5, 0.11), 0.2);
	
	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 1e6;
		i->k_s = 1e6;
		i->k_c = 1e6;
		i->density = 3000;
	}
}

void demTests::test3(){

	//Vaccuum
	vector4 unshocked(0.0, 0.0, 0.0, 0.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 0.2;
	domain.Lx = 4.0; domain.Ly = 4.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	//var.add_sphere(domain, Coordinates(0.5,0.5), 0.5);
	var.add_sphere(domain, Coordinates(1.5, 0.5), 0.5);
	var.add_sphere(domain, Coordinates(2.5, 0.5), 0.5);
	var.add_sphere(domain, Coordinates(3.5, 0.5), 0.5);

	var.gravity = 10.0;
	var.add_sphere(domain, Coordinates(2.0, 2.0), 0.5);
	var.particles.back().set_velocity(vector2(0.0, -2.0), 0.0);
	//var.add_sphere(domain, Coordinates(1.2, 0.3), 0.2);
	//var.add_sphere(domain, Coordinates(1.5, 0.11), 0.2);
	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 100;
		i->k_s = 100;
		i->k_c = 100;
		i->density = 0.625;
	}
}

void demTests::test4(){

	//Shock location
	//double x_s = 0.3;

	//------------------------------------------------------------
	//initial conditions for shock and unshocked fluid
	//vector4 shocked(1.3764, 0.394, 0.0, 1.5698); //density, velocity, pressure
	vector4 unshocked(0.0, 0.0, 0.0, 0.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 1.0;
	domain.Lx = 2.0; domain.Ly = 1.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------

	var.add_sphere(domain, Coordinates(0.7,0.2), 0.1);
	var.particles[0].set_velocity(vector2(0.0, 0.0), 0.0);

	var.add_sphere(domain, Coordinates(0.8, 0.3), 0.1);
	var.particles[1].set_velocity(vector2(0.0, 0.0), 0.0);

	var.add_sphere(domain, Coordinates(0.9, 0.2), 0.1);
	var.particles[2].set_velocity(vector2(0.0, 0.0), 0.0);

	var.add_sphere(domain, Coordinates(1.0, 0.3), 0.1);
	var.particles[3].set_velocity(vector2(0.0, 0.0), 0.0);

	var.add_sphere(domain, Coordinates(1.1, 0.2), 0.1);
	var.particles[4].set_velocity(vector2(0.0, 0.0), 0.0);

	var.add_sphere(domain, Coordinates(1.2, 0.3), 0.1);
	var.particles[5].set_velocity(vector2(0.0, 0.0), 0.0);

	var.add_sphere(domain, Coordinates(1.3, 0.2), 0.1);
	var.particles[6].set_velocity(vector2(0.0, 0.0), 0.0);

}

void demTests::test5(){
	//smooth acceleration (need to change subcycling to forced motion)
	//vector4 unshocked(2.6069, 0.6944, 0.0, 2.4583);
	vector4 unshocked(1.4, 0.0, 0.0, 1.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 12.0;
	domain.Lx = 4.0; domain.Ly = 2.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	var.add_sphere(domain, Coordinates(1.0, 1.0), 0.5);
	var.particles.back().set_velocity(vector2(0.0, 0.0), 0.0);
	var.particles.back().set_density(3.75);
}

void demTests::test6(){
	//2 particle shock

	vector4 shocked(2.6069, 0.6944, 0.0, 2.4583);
	vector4 unshocked(1.4, 0.0, 0.0, 1.0);
	//physical values
	//vector4 shocked(1.3333, 0.3535*sqrt(1e5), 0.0, 1.5*1e5);
	//vector4 unshocked(1.0, 0.0, 0.0, 1e5);

	initialL = shocked;
	initialR = unshocked;

	//domain.tstop = 0.5;
	domain.tstop = 1.5;
	domain.Lx = 4.0; domain.Ly = 4.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (i*domain.dx >= 0.4) interface(i, j) = true;
			else interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	//const double pi = boost::math::constants::pi<double>();
	//var.add_sphere(domain, Coordinates(1.0, 2.0), 0.5);
	//var.particles.back().set_density(1./(0.5*pi));

	Polygon poly;
	try{
		poly.create_from_file(domain,"vertices_zero.txt", vector2(1.5, 2.0));
		//poly.create_square(domain, 1.0, vector2(1.5, 2.0));
	}
	catch (const char c){
		std::cout << c << std::endl;
	}
	var.add_particle(poly, domain);
	

	//var.add_sphere(domain, Coordinates(2.22, 2.0), 0.2);
	//var.particles.back().set_velocity(vector2(0.0, 0.0), 0.0);
	//var.particles.back().set_density(1./(0.5*pi));
}

void demTests::test7(){
	//2 particle shock

	vector4 shocked(2.6069, 0.6944, 0.0, 2.4583);
	vector4 unshocked(1.4, 0.0, 0.0, 1.0);
	//physical values
	//vector4 shocked(1.3333, 0.3535*sqrt(1e5), 0.0, 1.5*1e5);
	//vector4 shocked(1.8621, 0.6944, 0.0, 2.4583*1e5);
	//vector4 unshocked(1.0, 0.0, 0.0, 1e5);

	initialL = shocked;
	initialR = unshocked;

	//domain.tstop = 0.5;
	domain.tstop = 2.0;
	domain.Lx = 4.0; domain.Ly = 4.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
			if (i*domain.dx >= 0.9) interface(i, j) = true;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	/*Polygon poly;
	try{
		poly.create_from_file(domain,"vertices_zero.txt", vector2(1.5, 2.0));
	}
	catch (const char c){
		std::cout << c << std::endl;
	}
	var.add_particle(poly, domain);
	*/
/*
	
	var.add_sphere(domain, Coordinates(1.5, 2.0), 0.5);
	var.particles.back().set_density(1./(0.5*pi));
	var.particles.back().set_mass(0.5);

	var.add_sphere(domain, Coordinates(2.6, 2.0), 0.5);
	var.particles.back().set_density(1./(0.5*pi));
	var.particles.back().set_mass(0.5);
*/
	const double pi = boost::math::constants::pi<double>();

	for (int i=1; i<4; i++){
		for (int j=0; j<2; j++){
			var.add_sphere(domain, Coordinates(0.98*i+0.5, 0.98*j+0.49), 0.5);
		}
	}

	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 100;
		i->k_s = 100;
		i->k_c = 100;
		i->density = 1./(0.5*pi);
	}

}

void demTests::test8(){
	//2 particle shock

	vector4 shocked(2.6069, 0.6944, 0.0, 2.4583);
	vector4 unshocked(1.4, 0.0, 0.0, 1.0);
	//physical values
	//vector4 shocked(1.3333, 0.3535*sqrt(1e5), 0.0, 1.5*1e5);
	//vector4 shocked(1.8621, 0.6944, 0.0, 2.4583*1e5);
	//vector4 unshocked(1.0, 0.0, 0.0, 1e5);

	initialL = shocked;
	initialR = unshocked;

	//domain.tstop = 0.5;
	domain.tstop = 2.0;
	domain.Lx = 4.0; domain.Ly = 4.0;
	domain.dx = domain.Lx/(domain.Nx-1);
	domain.dy = domain.Ly/(domain.Ny-1);
	domain.X.resize(domain.Nx, domain.Ny);

	var.fluid.state_function = StateFunctions::create(EOS_IG);
	var.fluid.state_function->y = 1.4;
	//------------------------------------------------------------
	//	Domain
	//------------------------------------------------------------
	//assigning x values
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			Coordinates point(i*domain.dx, j*domain.dy);
			domain.X(i, j) = point;
		}
	}
	interface.resize(domain.Nx, domain.Ny);
	//------------------------------------------------------------
	//	Fluid
	//------------------------------------------------------------
	//initial conditions of the fluid surrounding the levelset
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			interface(i, j) = false;
			if (i*domain.dx >= 0.9) interface(i, j) = true;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	/*Polygon poly;
	try{
		poly.create_from_file(domain,"vertices_zero.txt", vector2(1.5, 2.0));
	}
	catch (const char c){
		std::cout << c << std::endl;
	}
	var.add_particle(poly, domain);
	*/
/*
	
	var.add_sphere(domain, Coordinates(1.5, 2.0), 0.5);
	var.particles.back().set_density(1./(0.5*pi));
	var.particles.back().set_mass(0.5);

	var.add_sphere(domain, Coordinates(2.6, 2.0), 0.5);
	var.particles.back().set_density(1./(0.5*pi));
	var.particles.back().set_mass(0.5);
*/
	const double pi = boost::math::constants::pi<double>();

	for (int i=1; i<4; i++){
		for (int j=0; j<2; j++){
			var.add_sphere(domain, Coordinates(0.98*i+0.5, 0.98*j+0.49), 0.5);
		}
	}

	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 100;
		i->k_s = 100;
		i->k_c = 100;
		i->density = 1./(0.5*pi);
	}

}
//--------------------------------------------------------------------------------------------
//	Ghost Fluid Tests
//--------------------------------------------------------------------------------------------
/*void gfmTests::test1(){
	vector Left(1.0, 0.0, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.25;
	L = 1.0;

	var1.N = N;
	var1.dx = domain.L/domain.N;
	//var1.x0 = x0; //The GFM class has its own initial condition method
		//that deals directly with the gfm test struct
	//var1.tstop = tstop;
	var1.state_function = StateFunctions::create(EOS_IG);
	var1.state_function->y = 1.4;

	var2.N = N;
	var2.dx = domain.L/domain.N;
	//var2.x0 = x0;
	//var2.tstop = tstop;
	var2.state_function = StateFunctions::create(EOS_IG);
	var2.state_function->y = 1.4;
}

void gfmTests::test2(){ //Two rarefraction waves, vaccum test
	vector Left(1.0, -2.0, 0.4);
	vector Right(1.0, 2.0, 0.4);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.15;
	L = 1.0;
	yL = 1.4;
	yR = 1.4;
}

void gfmTests::test3(){ //strong shock wave of shock Mach number 198
	vector Left(1.0, 0.0, 1000.0);
	vector Right(1.0, 0.0, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.012;
	L = 1.0;
	yL = 1.4;
	yR = 1.4;
}

void gfmTests::test4(){ //Three strong discontinuities travelling to the right.
	vector Left(5.99924, 19.5975, 460.894);
	vector Right(5.99242, -6.19633, 46.0950);

	initialL = Left;
	initialR = Right;

	x0 = 0.4;
	tstop = 0.035;
	L = 1.0;
	yL = 1.4;
	yR = 1.4;
}

void gfmTests::test5(){ //slowly moving contact discontinuities
	vector Left(1.0, -19.59745, 1000.0);
	vector Right(1.0, -19.59745, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.8;
	tstop = 0.012;
	L = 1.0;
	yL = 1.4;
	yR = 1.4;
}

void gfmTests::test_example_1(){
	vector Left(2.0, 0.0, 9.8e5);
	vector Right(1.0, 0.0, 2.45e5);

	initialL = Left;
	initialR = Right;

	L = 4.0;
	x0 = 0.5;
	tstop = 0.0022;

	yL = 1.4;
	yR = 1.4;
}

void gfmTests::testA(){
	vector Left(1.0, 0.0, 1e5);
	vector Right(0.125, 0.0, 1e4);

	initialL = Left;
	initialR = Right;

	L = 1.0;
	x0 = 0.5;
	tstop = 0.0007;
	//tstop = 0;

	yL = 1.4;
	yR = 1.2;
}


void gfmTests::testB(){ //FedKiw 2002 test B
	number_of_materials = 3;
	vector Left(1.3333, 0.3535*sqrt(1e5), 1.5e5);
	vector Middle(1.0, 0.0, 1e5);
	vector Right(0.1379, 0.0, 1e5);

	initialL = Left;
	initialM1 = Middle;
	initialR = Right;

	L = 1.0;
	//N = 100;
	x0 = 0.05; //Right going shock between 5th ans 6th grid point (on a 100 grid)
	x1 = 0.5; //Material discontinuity between 50 and 51st point
	tstop = 0.0012;

	yL = 1.4;
	yR = 1.67;
	yM1 = 1.4;
}

void gfmTests::testC(){
	number_of_materials = 3;
	vector Left(1.3333, 0.3535*sqrt(1e5), 1.5e5);
	vector Middle(1.0, 0.0, 1e5);
	vector Right(3.1538, 0, 1e5);

	initialL = Left;
	initialM1 = Middle;
	initialR = Right;

	L = 1.0;
	x0 = 0.05; 
	x1 = 0.5; 
	tstop = 0.0017;

	yL = 1.4;
	yR = 1.249;
	yM1 = 1.4;
}

void gfmTests::testB_Wang(){
	number_of_materials = 4;
	vector Left(1.3333, 0.3535*sqrt(1e5), 1.5e5);
	vector MiddleLeft(1.0, 0.0, 1e5);
	vector MiddleRight(0.1379, 0, 1e5);
	vector Right(1.0, 0.0, 1e5);

	initialL = Left;
	initialM1 = MiddleLeft;
	initialM2 = MiddleRight;
	initialR = Right;

	L = 1.0;
	x0 = 0.05; 
	x1 = 0.4; 
	x2 = 0.6;
	tstop = 0.0014;


	yL = 1.4;
	yM1 = 1.4;
	yM2 = 1.67;
	yR = 1.4;
}

void gfmTests::testSG(){ //Water - Air shock tube
	stiffgas = true;
	vector Left(1000, 0.0, 1e9); //Water
	vector Right(50, 0.0, 1e5); //Air

	initialL = Left;
	initialR = Right;

	x0 = 0.7;
	tstop = 0.00023744;

	yL = 4.4;
	yR = 1.4;

	Pref1 = 6e8;
	Pref2 = 0.0; 
}

void gfmTests::testMach10(){
	//Wang test B with a mach 10 shock propagating right
	number_of_materials = 4;

	vector Left(5.92593, 6220.51, 4.665e7);
	vector MiddleLeft(1.0, 0.0, 1e5);
	vector MiddleRight(0.1379, 0, 1e5);
	vector Right(1.0, 0.0, 1e5);

	initialL = Left;
	initialM1 = MiddleLeft;
	initialM2 = MiddleRight;
	initialR = Right;

	L = 1.0;
	x0 = 0.05; 
	x1 = 0.4; 
	x2 = 0.6;

	tstop = 0.0001;

	yL = 1.4;
	yM1 = 1.4;
	yM2 = 1.67;
	yR = 1.4;
}

void gfmTests::testMach10_2(){ //testing it in Fedkiw's case

	number_of_materials = 3;
	vector Left(5.92593, 6220.51, 4.665e7);
	vector Middle(1.0, 0.0, 1e5);
	vector Right(0.1379, 0.0, 1e5);

	initialL = Left;
	initialM1 = Middle;
	initialR = Right;

	L = 1.0;
	//N = 100;
	x0 = 0.05; //Right going shock between 5th ans 6th grid point (on a 100 grid)
	x1 = 0.5; //Material discontinuity between 50 and 51st point
	tstop = 0.0002;

	yL = 1.4;
	yR = 1.67;
	yM1 = 1.4;
}

void gfmTests::testcase2(){
	number_of_materials = 2;
	yL = 5./3.;
	yR = 1.2;
	double yratio = (yL+1)/(yL-1);

	double d2 = (yratio*100+1)/(yratio+100)*0.82369077;
	vector Left(d2, 9.94949, 100.0);
	vector Right(1.0, 0.0, 1.0);

	initialL = Left;
	initialR = Right;

	L = 1.0;
	x0 = 0.2;
	tstop = 0.06;
}

//selecting the tests
void gfmTests::switch_resolution(){
	std::cout << "Resolution options" << std::endl
		<< "1. Low (100 Cells)" << std::endl
		<< "2. Mediumn (200 Cells)" << std::endl
		<< "3. High (400 Cells)" << std::endl
		<< "4. Very High (1000 Cells)" << std::endl
		<< "5. Exit" << std::endl;

	int a = Get_Switch();

	switch(a){
		case 0:
		case 1:
			N = 100;
			break;
		case 2:
			N = 200;
			break;
		case 3:
			N = 400;
			break;
		case 4:
			N = 1000;
			break;
		case 5:
		case 6:
		case 7:
		case 8:
		case 9:
			exit(0);
	}
}

void gfmTests::switch_test(){
	std::cout << "Test options" << std::endl << std::endl
		<< "Single Material Shock-Tube Tests" << std::endl
		<< "1. Test 1" << std::endl
		<< "2. Test 2" << std::endl
		<< "3. Test 3" << std::endl
		<< "Multimaterial Tests" << std::endl
		<< "4. Test B (Fedkiw 2002)" << std::endl
		<< "5. Test B (Wang 2004)" << std::endl
		<< "6. Mach 10 Test" << std::endl
		<< "7. Water-Air Test" << std::endl
		//<< "8. Test Problem 2" << std::endl
		<< "8. Exit" << std::endl;

	int a = Get_Switch();

	switch(a){
		case 0:
		case 1:
			test1();
			break;
		case 2:
			test2();
			break;
		case 3:
			test3();
			break;
		case 4:
			testB();
			break;
		case 5:
			testB_Wang();
			break;
		case 6:
			testMach10();
			break;
		case 7:
			testSG();
			break;
		case 8:
			//testcase2();
			//break;
		case 9:
			exit(0);
	}
}
*/







