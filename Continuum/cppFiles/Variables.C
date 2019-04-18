#include "../headerFiles/Variables.h"

//Polygon::~Polygon() {}

void Coordinates::scale(double xscale, double yscale){
	x = x * xscale;
	y = y * yscale;
}

void Coordinates::move(Coordinates point){
	x = x + point.x;
	y = y + point.y;
}

void Coordinates::display(){
	std::cout << '(' << x << ", " << y << ')' << '\t';
}

Coordinates Coordinates::operator +(const Coordinates &point){
	return Coordinates(x+point.x, y+point.y);
}

Coordinates Coordinates::operator -(const Coordinates &point){
	return Coordinates(x-point.x, y-point.y);
}

Coordinates Coordinates::operator *(const double &scalefactor){
	return Coordinates(x*scalefactor, y*scalefactor);
}

Coordinates Coordinates::operator /(const double &scalefactor){
	return Coordinates(x/scalefactor, y/scalefactor);
}



double Coordinates::length(){
	return sqrt(x*x + y*y);
}

void Polygon::create_square(Domain2D domain, double length, Coordinates center){
	n = 4;
	vertices.resize(n);
	edges.resize(n);

	//Generating a 1x1 square with vertices counterclockwise 
	//centered at the origin 
	vertices[0] = Vertex(0.5, 0.5);
	vertices[1] = Vertex(-0.5, 0.5);
	vertices[2] = Vertex(-0.5, -0.5);
	vertices[3] = Vertex(0.5, -0.5);

	for (int i=0; i<n; i++){
		edges[i].head = &vertices[i];
		if (i+1 < n) edges[i].tail = &vertices[i+1];
		else edges[i].tail = &vertices[0];
	}

	//displacement from origin = center(x, y)
	//scaling and moving the square
	for (int i=0; i<n; i++){
		vertices[i].point.scale(length, length);
		vertices[i].point.move(center);
	}

}

void Polygon::create(){

}

void Polygon::generate_surfacepoints(Domain2D domain){

	//Storing the index of points on the surface of polygon
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){

		}
	}
}

std::vector<Pos_Index> Bresenham::steep_pos(Domain2D domain, Coordinates P1, Coordinates P2){
	//stepping through in y (positive direction)
	double e = 0; //error
	double k = (P2.x - P1.x) / (P2.y - P1.y);

	//obtaining initial indices
	int j = ceil(P1.y/domain.dy); //j*dy must be greater than P1.y
	double y = j*domain.dy;

	//actual value of x at cell i 
	double x = k*(y - P1.y) + P1.x;
	int i = ceil(x/domain.dx); //j*dy > y
	if (i*domain.dx - x > 0.5*domain.dx) i = i-1;
	e = i*domain.dx - x;
	x = i*domain.dx;

	//indices storage
	std::vector<Pos_Index> indices;

	while (y <= P2.y){
		Pos_Index P(i, j);
		indices.push_back(P);

		if ((e + k*domain.dy) < 0.5*domain.dx){
			//point (i+1, j) is chosen
			e = e + k*domain.dy;
			j += 1;
		}

		else {
			e = e + k*domain.dy - domain.dx; 
			j += 1;
			i += 1;
		}

		x = i*domain.dx;
		y = j*domain.dy;
	}

	return indices;
}

std::vector<Pos_Index> Bresenham::gradual_pos(Domain2D domain, Coordinates P1, Coordinates P2){
	//stepping through in x
	double e = 0; //error
	double m = (P2.y - P1.y) / (P2.x - P1.x);

	//obtaining initial indices
	int i = ceil(P1.x/domain.dx); //i*dx must be greater than P1.x
	double x = i*domain.dx;

	//actual value of y at cell i 
	double y = m*(x - P1.x) + P1.y;
	int j = ceil(y/domain.dy); //j*dy > y
	if (j*domain.dy - y > 0.5*domain.dy) j = j-1;
	e = j*domain.dy - y;
	y = j*domain.dy;

	//indices storage
	std::vector<Pos_Index> indices;

	while (x <= P2.x){
		Pos_Index P(i, j);
		indices.push_back(P);

		if ((e + m*domain.dx) < 0.5*domain.dy){
			//point (i+1, j) is chosen
			e = e + m*domain.dx;
			i += 1;
		}

		else {
			e = e + m*domain.dx - domain.dy; 
			i += 1;
			j += 1;
		}

		x = i*domain.dx;
		y = j*domain.dy;
	}

	return indices;
}

std::vector<Pos_Index> Bresenham::steep_neg(Domain2D domain, Coordinates P1, Coordinates P2){
	//stepping through in y (negative direction)
	double e = 0; //error
	double k = (P2.x - P1.x) / (P2.y - P1.y);

	//obtaining initial indices
	int j = floor(P1.y/domain.dy); //j*dy must be less than P1.y
	double y = j*domain.dy;

	//actual value of x at cell i 
	double x = k*(y - P1.y) + P1.x;
	int i = ceil(x/domain.dx); //j*dy > y
	if (i*domain.dx - x > 0.5*domain.dx) i = i-1;
	e = i*domain.dx - x;
	x = i*domain.dx;

	//indices storage
	std::vector<Pos_Index> indices;

	while (y >= P2.y){
		Pos_Index P(i, j);
		indices.push_back(P);

		if ((e + k*domain.dy) > -0.5*domain.dx){
			//point (i+1, j) is chosen
			e = e + k*domain.dy;
			j -= 1;
		}

		else {
			e = e + k*domain.dy + domain.dx; 
			j -= 1;
			i += 1;
		}

		x = i*domain.dx;
		y = j*domain.dy;
	}

	return indices;
}

std::vector<Pos_Index> Bresenham::gradual_neg(Domain2D domain, Coordinates P1, Coordinates P2){
	//stepping through in x
	double e = 0; //error
	double m = (P2.y - P1.y) / (P2.x - P1.x); //m is negative

	//obtaining initial indices
	int i = ceil(P1.x/domain.dx); //i*dx must be greater than P1.x
	double x = i*domain.dx;

	//actual value of y at cell i 
	double y = m*(x - P1.x) + P1.y;
	int j = floor(y/domain.dy); //Initially assume the lower cell ,j*dy < y
	if (j*domain.dy - y < -0.5*domain.dy) j = j+1;
	e = j*domain.dy - y;
	y = j*domain.dy;

	//indices storage
	std::vector<Pos_Index> indices;

	while (x <= P2.x){
		Pos_Index P(i, j);
		indices.push_back(P);

		if ((e + m*domain.dx) > -0.5*domain.dy){
			//point (i+1, j) is chosen
			e = e + m*domain.dx;
			i += 1;
		}

		else {
			e = e + m*domain.dx + domain.dy; 
			i += 1;
			j -= 1;
		}

		x = i*domain.dx;
		y = j*domain.dy;
	}

	return indices;
}

std::vector<Pos_Index> Bresenham::line_algorithm(Domain2D domain, Coordinates P1, Coordinates P2){
	//P2 must always be to the right of P1
	//swap the two points if P2 is to the left of P1
	if (P2.x < P1.x){
		Coordinates tmp;
		tmp = P2;
		P2 = P1;
		P1 = tmp;
	}

	//calculate the slope
	double slope = (P2.y - P1.y) / (P2.x - P1.x);
	//case A, steep positive slope
	if (slope > 1){
		return steep_pos(domain, P1, P2);
	}
	//case B, gentle positive slope
	else if (slope > 0 && slope <= 1){
		return gradual_pos(domain, P1, P2);
	}
	//case C, steep negative slope
	else if (slope < -1){
		return steep_neg(domain, P1, P2);
	}
	//case D, gentle negative slope
	else {
		return gradual_neg(domain, P1, P2);
	}
}







