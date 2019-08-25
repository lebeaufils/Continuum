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
	//2 particles collision
	vector4 unshocked(0.0, 0.0, 0.0, 0.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 0.2;
	domain.Lx = 1.0; domain.Ly = 1.0;
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
	var.add_sphere(domain, Coordinates(0.3,0.5), 0.1);
	var.particles.back().set_velocity(vector2(0, 4), 0);
	var.add_sphere(domain, Coordinates(0.7,0.5), 0.1);
	var.particles.back().set_velocity(vector2(-2, 0), 0);
	/*
	Polygon poly;
	poly.create_from_file(domain,"vertices_zero.txt", vector2(1.5, 2.0));
	var.add_particle(poly, domain);
	var.particles.back().set_velocity(vector2(2, 0), 0);

	Polygon poly1;
	poly1.create_from_file(domain,"vertices_zero.txt", vector2(2.5, 2.0));
	var.add_particle(poly1, domain);
	var.particles.back().set_velocity(vector2(-2, 0), 0);
	*/
	

	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 200;//1e6;
		i->k_s = 200;//1e6;
		i->k_c = 200;//1e6;
		i->density = 0.625; //3000;
	}
}

void demTests::test2(){

	//vaccuum
	//12 particles with gravity
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
	var.gravity = 9.81;
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

	var.add_sphere(domain, Coordinates(0.6, 0.1), 0.1);
	var.add_sphere(domain, Coordinates(0.8, 0.1), 0.1);
	var.add_sphere(domain, Coordinates(1.0, 0.1), 0.1);
	var.add_sphere(domain, Coordinates(1.2, 0.1), 0.1);
	var.add_sphere(domain, Coordinates(1.4, 0.1), 0.1);

	var.add_sphere(domain, Coordinates(0.7, 0.32), 0.1);
	var.add_sphere(domain, Coordinates(0.9, 0.32), 0.1);
	var.add_sphere(domain, Coordinates(1.1, 0.32), 0.1);
	var.add_sphere(domain, Coordinates(1.3, 0.32), 0.1);

	var.add_sphere(domain, Coordinates(0.8, 0.55), 0.1);
	var.add_sphere(domain, Coordinates(1.0, 0.55), 0.1);
	var.add_sphere(domain, Coordinates(1.2, 0.55), 0.1);


	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 10.0;
		i->k_n = 1e6;
		i->k_s = 7e5;
		i->k_c = 1e6;
		i->density = 1000;
	}

}

void demTests::test3(){

	//flow
	//smooth acceleration (need to change subcycling to forced motion)
	vector4 unshocked(1.4, 0.0, 0.0, 1.0);

	initialL = unshocked;
	initialR = unshocked;

	domain.tstop = 8.0;
	domain.Lx = 3.0; domain.Ly = 2.5;
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
	var.add_sphere(domain, Coordinates(1.0, 1.25), 0.5);
	var.particles.back().set_density(0.625);
}

void demTests::test4(){
	//Flow around polygon

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
	for (int i=0; i<domain.Nx; i++){//
		for (int j=0; j<domain.Ny; j++){
			if (i*domain.dx >= 0.5) interface(i, j) = true;
			else interface(i, j) = false;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	const double pi = boost::math::constants::pi<double>();
	//var.add_sphere(domain, Coordinates(1.5, 2.0), 0.5);
	//

	Polygon poly;
	poly.create_from_file(domain,"vertices_zero.txt", vector2(1.5, 2.0));
	//poly.create_square(domain, 1.0, vector2(1.5, 2.0));
	var.add_particle(poly, domain);
	var.particles.back().set_density(1./(0.5*pi));
	

	//var.add_sphere(domain, Coordinates(2.22, 2.0), 0.2);
	//var.particles.back().set_velocity(vector2(0.0, 0.0), 0.0);
	//var.particles.back().set_density(1./(0.5*pi));
}

void demTests::test5(){
	//flow
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
			interface(i, j) = false;
			if (i*domain.dx >= 0.8) interface(i, j) = true;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------
	const double pi = boost::math::constants::pi<double>();
	/*
	var.add_sphere(domain, Coordinates(1.38, 2.0), 0.5);
	var.particles.back().set_density(1./(0.5*pi));
	var.particles.back().set_mass(0.5);

	var.add_sphere(domain, Coordinates(2.4, 2.0), 0.5);
	var.particles.back().set_density(1./(0.5*pi));
	var.particles.back().set_mass(0.5);
	*/

	Polygon poly;
	poly.create_from_file(domain,"vertices_zero.txt", vector2(1.38, 2.0));
	var.add_particle(poly, domain);

	Polygon poly1;
	poly1.create_from_file(domain,"vertices_zero.txt", vector2(2.4, 2.0));
	var.add_particle(poly1, domain);

	

	for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 100;
		i->k_s = 100;
		i->k_c = 100;
		i->density = 1./(0.5*pi);
	}

}

void demTests::test6(){
	//System of particles

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
			if (i*domain.dx >= 0.5) interface(i, j) = true;
		}
	}
	//------------------------------------------------------------
	//	Rigid body (level set definition)
	//------------------------------------------------------------

	const double pi = boost::math::constants::pi<double>();

	var.add_sphere(domain, Coordinates(1.0, 1.0), 0.2);
	var.add_sphere(domain, Coordinates(0.8, 0.6), 0.15);

	//var.add_sphere(domain, Coordinates(1.2, 1.5), 0.3); //polygon

	var.add_sphere(domain, Coordinates(1.4, 1.2), 0.18);

	//var.add_sphere(domain, Coordinates(1.9, 1.5), 0.2);
	//var.add_sphere(domain, Coordinates(2.3, 1.1), 0.3); //polygon
	var.add_sphere(domain, Coordinates(1.2, 0.4), 0.25);

	//var.add_sphere(domain, Coordinates(1.7, 0.7), 0.3); //polygon
	var.add_sphere(domain, Coordinates(2.2, 0.4), 0.2);
	var.add_sphere(domain, Coordinates(2.6, 1.0), 0.3);
	

	Polygon poly;
	poly.create_from_file(domain,"vertices_zero.txt", vector2(1.0, 1.5));
	var.add_particle(poly, domain);
	Polygon poly1;
	poly1.create_from_file(domain,"vertices_zero1.txt", vector2(2.0, 1.2));
	var.add_particle(poly1, domain);
	Polygon poly2;
	poly2.create_from_file(domain,"vertices_zero2.txt", vector2(1.7, 0.7));
	var.add_particle(poly2, domain);


for (std::vector<Particle>::iterator i = var.particles.begin(); i!=var.particles.end(); i++){
		i->damping_coefficient = 0.0;
		i->k_n = 100;
		i->k_s = 100;
		i->k_c = 100;
		i->density = 5./(0.5*pi);//1./(0.5*pi);
	}

}





