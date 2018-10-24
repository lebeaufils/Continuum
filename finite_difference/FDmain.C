#include "FDschemes.H"

int main(void){

	LxF linadv(100, 1.0, 1.0, 0.9); //N, L, a, c
	linadv.initial_conditions_square();
	linadv.plotname(999);

	//loop conditions
	int count = 0;
	double t = 0;
	double tStop = 0.2;

	do{
		if ((t+linadv.dt) > tStop){
			linadv.dt = tStop - t;
		}

		linadv.lax_friedrichs();
		linadv.boundary_conditions();

		t += linadv.dt;
		count += 1;
		if (count%10== 0){
			linadv.plotname(t);
		}
	}while (t<tStop);

	std::cout << count << std::endl;
}

