#include "eulerTests.h"

void eulerTests::test1(){ //sod's tube test
	vector Left(1.0, 0.75, 1.0); //density, velocity, pressure
	vector Right(0.125, 0.0, 0.1);

	initialL = Left;
	initialR = Right;

	x0 = 0.3;
	tstop = 0.2;
}

void eulerTests::test2(){ //Two rarefraction waves, vaccum test
	vector Left(1.0, -2.0, 0.4);
	vector Right(1.0, 2.0, 0.4);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.15;
}

void eulerTests::test3(){ //strong shock wave of shock Mach number 198
	vector Left(1.0, 0.0, 1000.0);
	vector Right(1.0, 0.0, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.5;
	tstop = 0.012;
}

void eulerTests::test4(){ //Three strong discontinuities travelling to the right.
	vector Left(5.99924, 19.5975, 460.894);
	vector Right(5.99242, -6.19633, 46.0950);

	initialL = Left;
	initialR = Right;

	x0 = 0.4;
	tstop = 0.035;
}

void eulerTests::test5(){ //slowly moving contact discontinuities
	vector Left(1.0, -19.59745, 1000.0);
	vector Right(1.0, -19.59745, 0.01);

	initialL = Left;
	initialR = Right;

	x0 = 0.8;
	tstop = 0.012;
}

/*--------------------------------------------------------------------------------------------
	Ghost Fluid Tests
--------------------------------------------------------------------------------------------*/
void gfmTests::set_number_of_cells(int Ncells){
	N = Ncells;
}

void gfmTests::test_example_1(){
	vector Left(2.0, 0.0, 9.8e5);
	vector Right(1.0, 0.0, 2.45e5);

	initialL = Left;
	initialR = Right;

	//L = 4.0;
	//N = 100;
	double dx = L/N;
	x0 = dx*(N+1)/2.;
	tstop = 0.0022;

	y1 = 1.4;
	y2 = 1.4;
}

void gfmTests::testA(){
	vector Left(1.0, 0.0, 1e5);
	vector Right(0.125, 0.0, 1e4);

	initialL = Left;
	initialR = Right;

	//L = 1.0;
	//N = 100;
	double dx = L/N;
	x0 = dx*(N+1)/2.;
	tstop = 0.0007;

	y1 = 1.4;
	y2 = 1.2;
}


void gfmTests::testB(){

}

void gfmTests::testC(){

}

void gfmTests::testD(){

}



