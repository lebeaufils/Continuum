#include "Functions.h"

int Functions::get_switch(){
	int switch_value;
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

int Functions::switch_resolution(){
	std::cout << "Resolution options" << std::endl
		<< "1. Low (100 Cells)" << std::endl
		<< "2. Mediumn (200 Cells)" << std::endl
		<< "3. High (400 Cells)" << std::endl
		<< "4. Exit" << std::endl;

	int a = get_switch();

	switch(a){
		case 0:
		case 1:
			return 100;
			break;
		case 2:
			return 200;
			break;
		case 3:
			return 400;
			break;
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
		case 9:
			exit(0);
	}
}

void Functions::get_double(double &x){
	std::cin >> x;
	//ignoring anything that cant be placed into the floating point number
	std::cin.ignore(32767, '\n'); 

	if (std::cin.fail()){
		std::cin.clear();
		std::cin.ignore(32767, '\n');
		std::cout << "Error, please enter a real number" << std::endl;
		get_double(x);
	}
}

void Functions::get_int(int &x){
	std::cin >> x;
	//ignoring anything that cant be placed into the floating point number
	std::cin.ignore(32767, '\n'); 

	if (std::cin.fail()){
		std::cin.clear();
		std::cin.ignore(32767, '\n');
		std::cout << "Error, please enter an integer" << std::endl;
		get_int(x);
	}
}