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

void standardTests::switch_resolution(){
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

//-----------------------------------------

void eulerTests::test1(){ //sod's tube test
	vector Left(1.0, 0.75, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.3;
	tstop = 0.2;

	var.N = N;
	var.dx = L/N;
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test1_stationary(){
	vector Left(1.0, 0.0, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.25;

	var.N = N;
	var.dx = L/N;
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test2(){ //Two rarefraction waves, vaccum test
	vector Left(1.0, -2.0, 0.4);
	vector Right(1.0, 2.0, 0.4);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.15;

	var.N = N;
	var.dx = L/N;
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test3(){ //strong shock wave of shock Mach number 198
	vector Left(1.0, 0.0, 1000.0);
	vector Right(1.0, 0.0, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.012;

	var.N = N;
	var.dx = L/N;
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test4(){ //Three strong discontinuities travelling to the right.
	vector Left(5.99924, 19.5975, 460.894);
	vector Right(5.99242, -6.19633, 46.0950);

	initialL = Left;
	initialR = Right;

	x0 = 0.4;
	tstop = 0.035;

	var.N = N;
	var.dx = L/N;
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;
}

void eulerTests::test5(){ //slowly moving contact discontinuities
	vector Left(1.0, -19.59745, 1000.0);
	vector Right(1.0, -19.59745, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.8;
	tstop = 0.012;

	var.N = N;
	var.dx = L/N;
	var.x0 = x0;
	var.tstop = tstop;
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
	tstop = 0.25;

	var.Nx = N;
	var.Ny = Ny;
	var.dx = L/(N-1);
	var.dy = Ly/(Ny-1);
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;

	for (int i=0; i<N; i++){
		x = i*var.dx;
		for (int j=0; j<Ny; j++){
			y = j*var.dy;
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
	tstop = 0.25;

	var.Nx = N;
	var.Ny = Ny;
	var.dx = L/(N-1);
	var.dy = Ly/(Ny-1);
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;

	for (int i=0; i<N; i++){
		x = i*var.dx;
		for (int j=0; j<Ny; j++){
			y = j*var.dy;
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
	tstop = 0.25;

	L = Ly = sqrt(0.5);

	var.Nx = N;
	var.Ny = Ny;
	var.dx = L/(N-1);
	var.dy = Ly/(Ny-1);
	var.x0 = x0;
	var.tstop = tstop;
	var.state_function = StateFunctions::create(EOS_IG);
	var.state_function->y = 1.4;

	//setting the list of interfacial points
	double x = 0;
	double y = 0;
	for (int i=0; i<N; i++){
		x = i*var.dx;
		for (int j=0; j<Ny; j++){
			y = j*var.dy;
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

//--------------------------------------------------------------------------------------------
//	Ghost Fluid Tests
//--------------------------------------------------------------------------------------------
/*
void gfmTests::set_number_of_cells(int Ncells){
	N = Ncells;
}

void gfmTests::test1(){
	vector Left(1.0, 0.0, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.25;
	L = 1.0;
	yL = 1.4;
	yR = 1.4;
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






