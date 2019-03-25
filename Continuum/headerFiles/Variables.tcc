//Template definitions

template <typename T>
T Euler2D::get_row(T A, int i){

	if (A.cols() == 0){
		throw "Array is empty.";
	}

	T B(A.cols, 1);
	for (int j=0; j<A.cols(); j++){
		B(j) = A(i, j);
	}
	return B;
}

template <typename T>
T Euler2D::get_column(T A, int i){

	if (A.rows() == 0){
		throw "Array is empty.";
	}
	
	T B(1, A.rows());
	for (int j=0; j<A.rows(); j++){
		B(j) = A(j, i);
	}
	return B;
}

template <typename T>
void Euler2D::display(T A){

	for (int i=0; i<A.rows(); i++){
		for (int j=(A.cols()-1); j>=0; j--){
			std::cout << A(i, j).transpose() << '\t';
		}
		std::cout << std::endl;
	}
}

template <typename T>
T Euler2D::swap_xy(T A){
	T B;
	B(0) = A(0);
	B(1) = A(3);
	B(2) = A(2);
	B(3) = A(1);
	return B;
}


