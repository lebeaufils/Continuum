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

void Domain2D::display_grid(){
	{
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				std::cout << X(i, j).x << ", " << X(i, j).y << '\t'; 
			}
			std::cout << std::endl;
		}
	}
}

int Polygon::orientation(Coordinates p0, Coordinates a, Coordinates b){
	//The sign of the cross product tells us if the points are oriented in an 
		// anti-clockwise orientation (negative)
		// clockwise orientation (positive)

	//Compute the cross product of vectors from the lowest point (determinant)
	double det = (a.x - p0.x) * (b.y - p0.y) - (b.x - p0.x) * (a.y - p0.y);
	if (det > 0)
		return -1;
    if (det < 0)
		return 1;
	return 0;
}

void Polygon::convex_hull(std::vector<Coordinates> &points){
	//O(nlogn)
	
	Coordinates p0 = points[0];
	int min_index = 0;

	//swap function
	auto swp = [](Coordinates &a, Coordinates &b){
		Coordinates tmp = a;
		a = b;
		b = tmp;
	};

	// l1-norm from p0
	auto dist = [p0](const Coordinates &a){
		return (a.x - p0.x)*(a.x - p0.x) + (a.y - p0.y)*(a.y - p0.y);
	};

	//find the lowest point
	for (int i=0; i < static_cast<int>(points.size()); i++){
		if (points[i].y < p0.y){
			p0 = points[i];
			min_index = i;
		}
		else if (points[i].y == p0.y){
			if (points[i].x < p0.x){
				p0 = points[i];
				min_index = i;
			}
		}
	}

	swp(points[0], points[min_index]);

	auto compare = [p0, dist](Coordinates a, Coordinates b){
	    int o = orientation(p0, a, b);
	    //if the points are colinear, priority is given to the nearer point
	    //this facilitates removal of colinear points later
	    if (o == 0)
	        return dist(a) <= dist(b);
	    //else, return true if a is to the left of b (anti clockwise)
	    return (o == -1);
	};

	std::sort(points.begin()+1,points.end(), compare);

	//now that a sorted list of points is obtained based on angle to the lowest point p0,
	//3 points are compared at a time, prev, current and next, to determine if the current point 
	//is part of the convex hull.

	//First, colinear points within the convex hull are removed.
	int m = 1;
	int n = static_cast<int>(points.size());
	for (int i=1; i<n; i++){
		//skipping over any points that have a colinear orientation
		while (orientation(p0, points[i], points[i+1]) == 0 && i<(n-1)){
			i++; //ignore all points that have have the same angle
			//the last point does not need to be compared,
			//and the comparison is halted at the second last point such that i+1 remains valid
		}

		//copying any points that do not have the same angle
		points[m] = points[i];
		m++; //size of the "new" array of points
	}

	std::ofstream outfile;
	outfile.open("points.txt");

	for (int i=0; i<m; i++){
		points[i].display();
		outfile << points[i].x << '\t' << points[i].y << std::endl;
	}
	std::cout << std::endl;
	outfile.close();

	std::cout << "m = " << m << std::endl;
	if (m<3) {
		throw "Convex hull does not exist for this set of points";
	}

	std::vector<Coordinates> chull;
	//pushing the first three points
	chull.push_back(points[0]);
	chull.push_back(points[1]);
	chull.push_back(points[2]);

	//searching the rest of the points
	//the top two points in the stack are compared with points in the list
	for (int i=3; i<m; i++){
		while (static_cast<int>(chull.size()) >= 2 && orientation(chull.end()[-2], chull.end()[-1], points[i]) != -1){
			std::cout << "Rejected Coordinates = ";
			chull.back().display(); std::cout << std::endl;
			chull.pop_back(); //remove the top value if it results in an anti-clockwise orientation
			//continue comparing the remaining points on the stack,
			//removing points until a clockwise orientation is obtained.
		}
		//once a left turn is established, push the next point 
		chull.push_back(points[i]);
	}

	vertices.clear();
	for (int i=0; i < static_cast<int>(chull.size()); i++){
		vertices.push_back(chull[i]);
	}
	

	//Polygon classification of the convex hull
	this->n = static_cast<int>(vertices.size());
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

std::vector<Coordinates> Polygon::random_points(double min, double max, int N){
	std::random_device r;
	std::seed_seq seed{ r(), r() };
	std::mt19937 rng(seed); //Mersenne Twister random engine  
	std::uniform_real_distribution<double> dist(min, max);

	std::vector<Coordinates> points;

	for(int i=0; i<N; i++){
		Coordinates p(dist(rng), dist(rng));
		points.push_back(p);
		//p.display();
	}

	return points;
}

void Polygon::output(Domain2D domain){
	std::ofstream outfile;
	std::ofstream outfile_2;

	outfile.open("polygon.txt");
	outfile_2.open("edgepoints.txt");

	for (int i=0; i<n; i++){
		outfile << vertices[i].x << '\t' << vertices[i].y << std::endl;
	}

	for (int i=0; i<static_cast<int>(surfacepoints.size()); i++){
		outfile_2 << domain.X(surfacepoints[i].i, surfacepoints[i].j).x <<
		'\t' << domain.X(surfacepoints[i].i, surfacepoints[i].j).y << std::endl;
	}

	/*boost::ptr_map<std::pair<int, int>, Edge>::iterator it = edges.begin();
	while (it != edges.end()){
		it->second->head.display();
		std::cout << '\t' << '\t';
		it->second->tail.display();
		std::cout << std::endl;
		it++;
	}*/

	outfile.close();
	outfile_2.close();
}

void Polygon::create_square(Domain2D domain, double length, Coordinates center){
	n = 4;

	//Generating a 1x1 square with vertices counterclockwise 
	//centered at the origin 
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
	//std::cout << edges.size() << std::endl;
	generate_surfacepoints(domain);
	output(domain);
}

void Polygon::create(Domain2D domain, double size, int K){
	//generates the convex hull of a random set of points
	double halflength_x = domain.Lx/2.;
	double halflength_y = domain.Ly/2.;
	double halflength = fmin(halflength_x, halflength_y);
	//creating a "uniformly" sized polygon
	//i.e. y and x have the same distribution
	if (size >= halflength){
		throw "Size of polygon exceeds domain size";
	}
	std::vector<Coordinates> points = random_points(halflength - size, halflength + size, K);
	try {
		convex_hull(points);
	}
	catch (const char c){
		std::cout << c << std::endl;
	}
	//std::cout << n << std::endl;
	generate_edges(vertices);
	//std::cout << edges.size() << std::endl;
	generate_surfacepoints(domain);
	output(domain);
}

int Polygon::point_in_polygon(Coordinates p){
	//PNPOLY Algorithm from Copyright (c) 1970-2003, Wm. Randolph Franklin
	//https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html#The%20Method
	
	//Bounding box
	double xmin = 1e6;
	double xmax = 0;
	double ymin = 1e6;
	double ymax = 0;

	for (std::vector<Vertex>::iterator it = vertices.begin(); it < vertices.end(); it++){
		if (it->x > xmax) xmax = it->x;
		if (it->x < xmin) xmin = it->x;
		if (it->y > ymax) ymax = it->y;
		if (it->y < ymin) ymin = it->y;
	}

	//Out of polygon boundaries
	if (p.x < xmin || p.x > xmax || p.y < ymin || p.y > ymax) {
    	return 0;
	}

	//Crossings test using the jordan curve
	int c = 0;
	//if ( ((verty[i]>testy) != (verty[j]>testy))
	//The first equality checks if the test point is within the y range of the edge
	//The second equality identifies if the point is to the left or right of the edgeline
	//if the line equation (point and slope) is satisfied, point p lies on the line
	//hence if p.x is less than this value, it lies to the left.
	//
	//effectively, an imaginary horizontal line is projected towards the right
	//if p lies to the left of the line, this imaginary ray will cross the line
	//the line equation is given by x - x0 = (x1 - x0)/(y1 - y0) * (y - y0)
	///////
	Edge* current_edge = edges.begin()->second;
	do{
		if ( ((current_edge->tail.y > p.y) != (current_edge->head.y > p.y)) &&
			(p.x < (current_edge->tail.x - current_edge->head.x) * (p.y - current_edge->head.y) / (current_edge->tail.y - current_edge->head.y) + current_edge->head.x) ){
			c = !c; //flips if it crosses an edge
		}
		current_edge = current_edge->next;
	}while(current_edge != edges.begin()->second);
	return c;
}

//
//Rastor representstion of lines
//
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

void LevelSet::display_grid(){
	for (int i=1; i<phi.rows()-1; i++){
		for (int j=1; j<phi.cols()-1; j++){
			std::cout << phi(i, j) << '\t';
		}
		std::cout << std::endl;
	}
}

void RB_2D::add_levelset(){
	LevelSet ls;
	levelsets.push_back(ls);
}







