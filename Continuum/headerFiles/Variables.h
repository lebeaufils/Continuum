#ifndef VARIABLES_H_
#define VARIABLES_H_

#include <Eigen/Dense>
#include <boost/ptr_container/ptr_map.hpp>
#include <vector>
#include <array>
//#include <map>
//#include <utility>
#include <cstdlib>
#include <stack> 
#include <random>
#include <fstream>
#include "EOS.h"


typedef Eigen::Vector3d vector;
typedef Eigen::Vector4d vector4;
typedef Eigen::MatrixXd matrix;
typedef Eigen::Matrix<vector4, Eigen::Dynamic, Eigen::Dynamic> vecarray;

//-----------------------------------------------------------------------------------------------
//	Data storage
//-----------------------------------------------------------------------------------------------
struct Coordinates
{
	double x;
	double y;

	Coordinates() : x(0), y(0) {}
	Coordinates(double x, double y) : x(x), y(y) {}
	Coordinates(const Coordinates &obj) {
   		x = obj.x;
   		y = obj.y;
	}

	void scale(double, double);
	void move(Coordinates);

	//overloaded operators
	Coordinates operator +(const Coordinates&);
	Coordinates operator -(const Coordinates&);
	Coordinates operator *(const double&);
	Coordinates operator /(const double&);

	double length();
	void display();
};

struct Pos_Index{
	int i; //x index
	int j; //y index
	//position indices for reference to a global array eg X(i, j)

	Pos_Index() : i(0), j(0) {}
	Pos_Index(int i, int j) : i(i), j(j) {}
};

struct Domain1D
{
	int N;
	double L;
	double dt;
	double dx;
	double tstop;

	matrix X;

	Domain1D() : N(0), L(0), dt(0), dx(0), tstop(0), X(0, 0) {}
	Domain1D(int N) : N(N), L(1.0), dt(0), dx(0), tstop(0), X(0, 0) {}
};

struct Domain2D
{
	int Nx;
	int Ny;
	double Lx;
	double Ly;
	double dt;
	double dx;
	double dy;
	double tstop;

	Eigen::Array<Coordinates, Eigen::Dynamic, Eigen::Dynamic> X; //coordinate pairs (x,y)

	Domain2D() : Nx(0), Ny(0), Lx(0), Ly(0), dt(0), dx(0), tstop(0), X(0, 0) {}
	Domain2D(int N) : Nx(N), Ny(N), Lx(1.0), Ly(1.0), dt(0), dx(0), tstop(0), X(N, N) {}
	Domain2D(int Nx, int Ny) : Nx(Nx), Ny(Ny), Lx(1.0), Ly(1.0), dt(0), dx(0), tstop(0), X(Nx, Ny) {}

	void display_grid(){
		for (int j=0; j<Ny; j++){
			for (int i=0; i<Nx; i++){
				std::cout << X(i, j).x << ", " << X(i, j).y << '\t'; 
			}
			std::cout << std::endl;
		}
	}
};

struct Euler1D
{
	//domain parameters
	//int N;
	//double dt;
	//double dx;
	//double x0;
	//double tstop;

	//Variable Matrix
	//matrix X; //Domain
	matrix U; //conserved variables
	matrix F; //flux

	std::shared_ptr<StateFunctions> state_function;

	//Euler1D() : N(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), state_function(NULL) {}
	Euler1D() : U(0, 0), F(0, 0), state_function(NULL) {}
	~Euler1D() {}
};


struct Euler2D
{
	//domain parameters
	//int Nx;
	//int Ny;
	//double dt;
	//double dx;
	//double dy;
	//double x0;
	//double tstop;

	//Variable Matrix
	//matrix X; //Domain
	vecarray U; //conserved variables
	vecarray F; //flux for the x derivative
	vecarray G; // "" y derivative

	std::shared_ptr<StateFunctions> state_function;

	//Euler2D() : Nx(0), Ny(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), G(0, 0), state_function(NULL) {}
	Euler2D() : U(0, 0), F(0, 0), G(0, 0), state_function(NULL) {}
	~Euler2D() {}

	//functions to return row i
	template <typename T>
	T get_row(T, int);
	template <typename T>
	T get_column(T, int);
	template<typename T>
	void display(T);
	template<typename T>
	T swap_xy(T); //swaps the order of velocity in x-y directions
};

struct LevelSet
{
	matrix phi;

	LevelSet() : phi(0, 0) {}
};

struct RB_2D
{
	//collection of level sets, one for each object
	//Eigen::Array<matrix, 1, Eigen::Dynamic> levelset_array;
	std::vector<LevelSet> levelsets; //list of n-levelset
	Euler2D fluid;
	Euler2D rigidbody;

	RB_2D() : levelsets(0), fluid(), rigidbody() {}
};

//Polygon storage

struct Vertex : public Coordinates
{
	Vertex() : Coordinates(0, 0) {}
	Vertex(Coordinates xy) : Coordinates(xy) {} //constructor to copy base class
	Vertex(double x, double y) : Coordinates(x, y) {}
	~Vertex() {}

	//generate random polygons in the future.
};

struct Edge
{
	Vertex head;
	Vertex tail;

	Edge* next;

	//coordinates of edges are stored in domain struct as indices of 
	//the global array X(i, j)

	Edge() : head(), tail(), next(NULL) {}
	~Edge() {
		next = NULL;
	}

};

struct Polygon
{	
	//for an n-sided polygon
	int n;

	std::vector<Vertex> vertices;
	//std::map <std::pair<int, int>, Edge*> edges; //map of vertex pairs to corresponding edge
	boost::ptr_map<std::pair<int, int>, Edge> edges;
		//the pointers here are owned by the polygon class
	std::vector<Pos_Index> surfacepoints; //list of indices of points on the polygon surface

	Polygon() : n(0), edges(), surfacepoints(0) {}
	Polygon(int n) : n(n), edges(), surfacepoints(0) {}
	~Polygon() {}

	void convex_hull(std::vector<Coordinates>&); //graham scan

	void generate_edges(std::vector<Vertex>);
	void generate_surfacepoints(Domain2D);
	void output(Domain2D);

	void create_square(Domain2D, double, Coordinates);
	void create(Domain2D, double,int);

	int point_in_polygon(Coordinates);

	static int orientation(Coordinates, Coordinates, Coordinates);
	static std::vector<Coordinates> random_points(double, double, int);

};

struct Bresenham{
	static std::vector<Pos_Index> steep_pos(Domain2D, Coordinates, Coordinates);
	static std::vector<Pos_Index> gradual_pos(Domain2D, Coordinates, Coordinates);
	static std::vector<Pos_Index> steep_neg(Domain2D, Coordinates, Coordinates);
	static std::vector<Pos_Index> gradual_neg(Domain2D, Coordinates, Coordinates);

	static std::vector<Pos_Index> line_algorithm(Domain2D, Coordinates, Coordinates);
	static std::vector<Pos_Index> line_algorithm(Domain2D, Edge*);
};


#include "Variables.tcc"




#endif /* VARIABLES_H_ */