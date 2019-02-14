#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>

typedef Eigen::Vector3d vector;
typedef Eigen::MatrixXd matrix;

struct Functions
{
	//Functions that deal with user input and switching between different EOS, Solvers or tests.
	int get_switch();
	void switch_resolution();
	void switch_EOS();
	void switch_solver();
	void user_input();
	void settings_file();
};

int Functions::get_switch(){
	do{
		std::cin >> switch_value;

		bool find;
		if (switch_value >0 && switch_value < 10) find = 1;
		else find = 0;

		if (find == 1){
			int a = switch_value;
			break;
		}
		else{ //error if find fails
			std::cout << "Invalid input." << std::endl;
		}

	}while(true);
	return switch_value;
}

void Functions::switch_resolution(){
	std::cout << "Resolution options" << std::endl
		<< "1. Low (100 Cells)" << std::endl
		<< "2. Mediumn (200 Cells)" << std::endl
		<< "3. High (400 Cells)" << std::endl
		<< "4. Exit" << std::endl;

	int a = get_switch();

	switch(a){
	case 0:
	case 1:
		break;
	case 2:
		break;
	case 3:
		break;
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
	case 9:
		exit(0);
}


#endif /* FUNCTIONS_H_ */