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

int Polygon::orientation(Coordinates p0, Coordinates p1, Coordinates p2){
	//slope_p0p1 - slope_p2p1 = const
	int a = (p1.y - p0.y)*(p2.x - p1.x) - (p2.y - p1.y)*(p1.x - p0.x); 

	/*if (a > 0) return 1;
	else if (a < 0) return -1; 
	else return 0;*/

	return (a > 0) - (a < 0);
}

void Polygon::generate_edges(std::vector<Vertex> vertices){
	vertices.push_back(vertices[0]);

	if (!edges.empty()) {
		throw "There is a pre-existing list of edges";
	}
	//assuming vertices are sorted
	//creating the firstt edge node
	std::pair<int, int> v01 = std::make_pair(0, 1);
	edges.insert(v01, new Edge());
		edges[v01].head = vertices[0];
		edges[v01].tail = vertices[1];
		Edge* prevnode = &edges[v01];
		
	//create n edges with their assigned vertices
	for (int i=1; i<n; i++){
		std::pair<int, int> vij = std::make_pair(i, i+1);
		edges.insert(vij, new Edge());
		edges[vij].head = vertices[i];
		edges[vij].tail = vertices[i+1];

		prevnode->next = &edges[vij];
		prevnode = &edges[vij];
	}
	//link the final node to the first
	prevnode->next = &edges[v01];	
}

void Polygon::generate_surfacepoints(Domain2D domain){
	//iterating through the map of edges

	boost::ptr_map<std::pair<int, int>, Edge>::iterator i=edges.begin();
	while (i != edges.end()){
		Edge *edge = i->second; //i->second is the Edge* value in the map

		std::vector<Pos_Index> edgepoints = Bresenham::line_algorithm(domain, edge);

		surfacepoints.insert(surfacepoints.end(),
			std::make_move_iterator(edgepoints.begin()),
			std::make_move_iterator(edgepoints.end()));
		i++;
	}
}

void Polygon::create_square(Domain2D domain, double length, Coordinates center){
	n = 4;

	//Generating a 1x1 square with vertices counterclockwise 
	//centered at the origin 
	std::vector<Vertex> vertices;
	vertices.resize(n);

	vertices[0] = Vertex(0.5, 0.5);
	vertices[1] = Vertex(-0.5, 0.5);
	vertices[2] = Vertex(-0.5, -0.5);
	vertices[3] = Vertex(0.5, -0.5);

	//displacement from origin = center(x, y)
	//scaling and moving the square
	for (int i=0; i<n; i++){
		vertices[i].scale(length, length);
		vertices[i].move(center);
	}

	generate_edges(vertices);
	std::cout << edges.size() << std::endl;
	//generate_surfacepoints(domain);
}

void Polygon::create(){
	//generates the convex hull of a random set of points
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

std::vector<Pos_Index> Bresenham::line_algorithm(Domain2D domain, Edge* edge){
	//P2 must always be to the right of P1
	//swap the two points if P2 is to the left of P1
	Coordinates P2 = edge->tail;
	Coordinates P1 = edge->head;
	//Object slicing occurs here

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







