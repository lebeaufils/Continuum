#include "FVschemes.H"

int main(void){

	euler test(100, 1.0, 1.0, 0.9); //N, L, a, c
	test.initial_conditions_step();
	test.plotname(999);

	//loop conditions
	/*int count = 0;
	double t = 0;
	double tStop = 0.5;
	do{
		if ((t+LA.dt) > tStop){
			LA.dt = tStop - t;
		}

		//define flux vector
		LA.evaluate_flux();
		LA.solver();
		LA.boundary_conditions();

		t += LA.dt;
		count += 1;
		if (count%10 == 0){
			LA.plotname(t);
		}
	}while (t<tStop);*/
}






