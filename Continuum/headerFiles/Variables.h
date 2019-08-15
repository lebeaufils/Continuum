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
typedef Eigen::Vector2d vector2;
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
	Coordinates(const vector2 &obj) {
   		x = obj(0);
   		y = obj(1);
	}

	void scale(double, double);
	void move(Coordinates);

	//overloaded operators
	Coordinates operator +(const Coordinates&);
	Coordinates operator -(const Coordinates&);
	Coordinates operator *(const double&);
	Coordinates operator /(const double&);

	double length() const;
	void display() const;
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
	int buffer;//number of cells for extended domain to allow for overlap with boundaries during contact
	double Lx;
	double Ly;
	double dt;
	double dx;
	double dy;
	double tstop;

	Eigen::Array<Coordinates, Eigen::Dynamic, Eigen::Dynamic> X; //coordinate pairs (x,y)

	Domain2D() : Nx(0), Ny(0), buffer(0), Lx(0), Ly(0), dt(0), dx(0), tstop(0), X(0, 0) {}
	Domain2D(int N) : Nx(N), Ny(N), buffer(1+floor(N/20)), Lx(1.0), Ly(1.0), dt(0), dx(0), tstop(0), X(N, N) {}
	Domain2D(int Nx, int Ny) : Nx(Nx), Ny(Ny), buffer(1+floor((Nx+Ny)/40)), Lx(1.0), Ly(1.0), dt(0), dx(0), tstop(0), X(Nx, Ny) {}

	void display_grid() const;
};

struct Euler1D
{
	matrix U; //conserved variables
	matrix F; //flux

	std::shared_ptr<StateFunctions> state_function;

	//Euler1D() : N(0), dt(0), dx(0), x0(0), tstop(0), X(0, 0), U(0, 0), F(0, 0), state_function(NULL) {}
	Euler1D() : U(0, 0), F(0, 0), state_function(NULL) {}
	~Euler1D() {}
};


struct Euler2D
{
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
	LevelSet(int x, int y) : phi(x, y) {}

	void display_grid();
};

struct Stationary_RB
{
	//collection of level sets, one for each object
	//Eigen::Array<matrix, 1, Eigen::Dynamic> levelset_array;
	std::vector<LevelSet> levelsets; //list of n-levelset
	Euler2D fluid;
	LevelSet combinedls;
	Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal;

	Stationary_RB() : levelsets(0), fluid(), combinedls(), normal(0, 0) {}

	void add_levelset();
	void zeroes(int Nx, int Ny); //resize and set zeroes
};

//Polygon storage

struct Vertex : public Coordinates
{
	Vertex() : Coordinates(0, 0) {}
	Vertex(const Coordinates& xy) : Coordinates(xy) {} //constructor to copy base class
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
	void create_from_file(Domain2D);
	void create_from_file(Domain2D, std::string, vector2);

	int point_in_polygon(Coordinates) const;

	//void translate(double, double);

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

//------Rotors------
//modified from http://marctenbosch.com/quaternions/
//------------------
struct Bivector2
{
	//in 2D, the bivector of two vectors is confined to a single plane in xy
	double b12;

	Bivector2(double b) : b12(b) {}

	static Bivector2 exterior_product(const vector2&, const vector2&);
};

struct Rotor2
{	
	//exponential form
	// R = ab = cos(theta) + a^b
	//magnitude of rotation
	double a; 
	//bivector
	double b12; //unit bi-vector. in 2D, only one plane e12 exists and there is
	//only one choice for the unit plane

	//constructors
	Rotor2() : a(1), b12(1) {}
	Rotor2(double a, double b) : a(a), b12(b) {}
	//Rotor2(double a, const Bivector2 &bv) : a(a), b12(bv.b12) {}
	//advanced initialisation
	//Given initial and final vectors
	Rotor2(const vector2&, const vector2&);
	//Given angle and normalised plane(axis)
	Rotor2(double); //plane is bounded to xy in 2D
	Rotor2(double, const Bivector2&);
	~Rotor2() {};

	//rotate a vector
	vector2 rotate(const vector2& v) const;

	//utility functions
	Rotor2 operator*(const Rotor2&) const; //multiply 2 rotors
	Rotor2 operator/(const Rotor2&) const; //divide 2 rotors
	double sqlength() const; //cheaper computation if sqrt is not needed
	double length() const;
	void normalise();
	Rotor2 reverse() const; //reverses the direction of the bivector -- complex conjugatte
	Rotor2 nrotor() const; //normalised rotor
	Rotor2 inverse() const; //reverses the rotation
	//static Rotor2 geometric_product(const Rotor2&, const Rotor2&); //geometric product of two rotors
	//Rotor2 rotate_rotor(const Rotor2&);

	//wrapper to rotate about point p with a given angular velocity and lapsed time
	static vector2 rotate_about(const vector2&, const vector2&, double, double);
		static vector2 rotate_about(const vector2&, const vector2&, double);
	static vector2 rotate_reverse(const vector2&, const vector2&, double, double);
	static vector2 rotate_reverse(const vector2&, const vector2&, double);
};




#include "Variables.tcc"




#endif /* VARIABLES_H_ */