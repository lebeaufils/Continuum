#include "eulerTests.h"

void eulerTests::test1(){ //sod's tube test
	vector L(1.0, 0.75, 1.0); //density, velocity, pressure
	vector R(0.125, 0.0, 0.1);

	initialL = L;
	initialR = R;

	x0 = 0.3;
	tstop = 0.2;
}

void eulerTests::test2(){ //Two rarefraction waves, vaccum test
	vector L(1.0, -2.0, 0.4);
	vector R(1.0, 2.0, 0.4);

	initialL = L;
	initialR = R;

	x0 = 0.5;
	tstop = 0.15;
}

void eulerTests::test3(){ //strong shock wave of shock Mach number 198
	vector L(1.0, 0.0, 1000.0);
	vector R(1.0, 0.0, 0.01);

	initialL = L;
	initialR = R;

	x0 = 0.5;
	tstop = 0.012;
}

void eulerTests::test4(){ //Three strong discontinuities travelling to the right.
	vector L(5.99924, 19.5975, 460.894);
	vector R(5.99242, -6.19633, 46.0950);

	initialL = L;
	initialR = R;

	x0 = 0.4;
	tstop = 0.035;
}

void eulerTests::test5(){ //slowly moving contact discontinuities
	vector L(1.0, -19.59745, 1000.0);
	vector R(1.0, -19.59745, 0.01);

	initialL = L;
	initialR = R;

	x0 = 0.8;
	tstop = 0.012;
}


