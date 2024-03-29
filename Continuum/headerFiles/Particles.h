#ifndef PARTICLES_H_
#define PARTICLES_H_

#if   __cplusplus < 201103L
#error This program requires a C++11 compiler. 
#endif

#include "LevelSet.h"
#include <list> //stl linked list
#include <unordered_map> //to be replaced by perfect hash function if time permits
//#include <boost/functional/hash.hpp>
#include <memory>
#include <algorithm>

//------------------------------------------------------------
//Particles 
//------------------------------------------------------------
//Particles are represented by a levelset function
//The 
struct Particle //NEEDS WORK
{
	//RB_2D should store particles rather than levelsets.
	//each particle can have its reference levelset and nodes
	//------------------------------------------
	// Material Properties
	//------------------------------------------
	double mass;
	double density;
	double damping_coefficient; //Dampens the oscillation to simulate dissipation of energy
	double miu; //interparticle friction coefficient
	double k_n; //normal contact stiffness
	double k_s; //shear contact stiffness
	double k_c;

	//------------------------------------------
	// Location information //NOT IMPLEMENTED?
	//------------------------------------------
	int label; //particle label and resolution are assigned during hash table generation and will remain constant from then

	double size; //size of AABB, length of box
	int resolution; //minimum resolution required
	std::vector<std::pair<int, int> > igrids; //Intersected grids
	//Reference 
	LevelSet ls; //reference levelset - domain wide
	vector2 centroid; //initial, fixed
	std::vector<vector2> ref_nodes;

	//Dynamic
	LevelSet dynamicls; //levelset that moves with the particle - bucket wide
	vector2 centre; //tracking the center of dynamicls
	std::vector<vector2> nodes;

	//------------------------------------------
	// Kinematics
	//------------------------------------------
	vector2 vc; //translational velocity
	double w; //angular velocity
	vector2 s; //total displacement
	double theta; //total rotation

	//Forces
	vector2 force; //Linear force
	//One spring for EACH node!!
	std::vector<std::unordered_map<int, vector2> > springs; //Tangential spring that accumulates strain for each collision //analogous to the rows of close_tracker
	std::vector<std::unordered_map<int, vector2> > wall_springs; //Tangential spring that accumulates strain for wall collisions
	std::vector<std::unordered_map<int, double> > overlap;
	std::vector<std::unordered_map<int, double> > wall_overlap;
	double torque; //Accumuluator for torque

	//Tracking collisions? could be done with a global collision list?
	//bool in_collision; 

	Particle() : mass(1.0), density(3000), damping_coefficient(0.2), miu(0.26), k_n(1e5), k_s(1e5), k_c(1e5), label(0), size(0), resolution(0), igrids(0), ls(), centroid(0, 0), ref_nodes(0), dynamicls(), centre(0, 0), nodes(0), 
	vc(0, 0), w(0), s(0, 0), theta(0), force(0, 0), springs(0), wall_springs(), overlap(), wall_overlap(), torque(0) {}
	Particle(const Domain2D&, const Coordinates&, double);
	Particle(const Polygon&, const Domain2D&);
	//Particle(const Particle&) //copy constructor
	Particle(const Particle& gr);
	~Particle() {};

	//void initialise(const Polygon&, const Domain2D&);
	//Particle(double d); //no real need to specify density for rigid bodies
	double diameter(); //estimated particle diameter, equidiv?;

	void set_mass(double);
	void set_density(double);
	void set_velocity(const vector2&, double);
	LevelSet motion(const Domain2D&, const vector2&, double);

	//-----------------------------------------------------
	//Inertial properties
	//-----------------------------------------------------
	static double compute_mass(const Particle&, const Domain2D&);
	static vector2 center_of_mass(const LevelSet&, const Particle&, const Domain2D&);
	static double moment_of_inertia(const LevelSet&, const Particle&, const Domain2D&);
	static vector2 velocity(const Coordinates&, const Particle&); //vb

	static vector2 cross(double, const vector2&);
	static double cross(const vector2&, const vector2&);
	static LevelSet merge(const std::vector<Particle>&, const Domain2D&);
	static vector2 normal_sum(const LevelSet&, const Domain2D&);
	static double length_sum(const LevelSet&, const Domain2D&);
	//-----------------------------------------------------
	////Size of bounding box
	//-----------------------------------------------------
	//static double AABB_extent(const Polygon&); //using the vertices, for problems where polygons are generated by vertices
	static double AABB_extent(const Particle&); //using the level set function
};

//-----------------------------------------------------
// Hierarchical hash table
//-----------------------------------------------------
//Each divided cell in the domain is stored as a hash bucket, 
//with its bottom left corner bring the bucket's key.
//The value of this entry is a linked list containing the labels of all the
//particles that intersect the cell.
struct HashEntry
{
	std::pair<int, int> key;
	std::list<int> labels;

	//constructors
	HashEntry(std::pair<int, int> key_input) : key(key_input), labels(0) {}
	HashEntry(std::pair<int, int> key_input, int label);
	~HashEntry() {}

	//Accessing and modifying particle list
	void add_label(int);
	void remove_label(int);
	std::list<int> return_all_labels(int); //returns all labels except current
};

struct HashPair
{//hash function mapping grid location to an index 
//using the elegant paring method by Matthew Szudzik
	//by overloading (), we allow the class instance to be used as a function
	template <class a, class b>
	std::size_t operator()(const std::pair<a, b>& pair) const {
		auto hash = (pair.first < pair.second) ? pair.second * pair.second + pair.first : pair.first * pair.first + pair.first + pair.second;
    	return hash;
    }
	//equality checks are based on std::pair
};

struct HashTable
{
	std::unordered_multimap<std::pair<int, int>, int, HashPair> map;

	const int resolution; //resolution number is its position in the array at HHT
	const double resolution_size; //tile size 

	HashTable(int res, double res_size);
	~HashTable() {} 

	//-------------------------
	//Utility operations
	//-------------------------
	//individual particle insert and delete
	//void insert(const std::pair<int, int>&, int); //inserts a particle label into the key's list
	//void remove(const std::pair<int, int>&, int); //finds and remove a particle from the key's list
	//bool search(const std::pair<int, int>&, int); //returns all particles in the grid that are not the target particle
};

struct HierarchicalHashTable //HHT
{
	std::vector<HashTable> tables;
	std::vector<std::vector<int> > close_tracker; //2D array checking tracking closeness of each particle to the other, size nsquare
	std::vector<std::vector<int> > wall_tracker; //0-left, 1-right, 2-bottom, 3-top

	HierarchicalHashTable() : tables(), close_tracker(), wall_tracker() {}
	HierarchicalHashTable(int n);
	~HierarchicalHashTable() {}

	void create_close_tracker(int); //generates the close and wall trackers based on number of particles
	//generating table and adding layers to the hierarchy
	void add_table(double); //takes resolution size as input
	void generate_tables(const Domain2D&, std::vector<Particle>&);
	int compute_resolution(double);
	std::vector<std::pair<int, int> > find_intersecting_cells(double, const vector2&, int);
	void close_to_wall(const Domain2D&, const Particle&, int);
	void add_particle_all(const Domain2D&, const std::vector<Particle>&, const Particle&, int);

	//modifying  particle labels within the tables
	void move_particle(const Domain2D&, const std::vector<Particle>&, Particle&);
	void remove_particle(const std::vector<Particle>&, const Particle&, int, const std::vector<std::pair<int, int> >&);
	void add_particle(const std::vector<Particle>&, const Particle&, int, const std::vector<std::pair<int, int> >&);
	//When adding particle X, the bounding box is checked for intersection with the grid at resolution of X
};

struct Moving_RB
{
	//Rigid body system of moving particles
	//collection of "particles" which contain tHier own levelset and discretised nodes
	std::vector<Particle> particles; //A particle's label is its position in this vector, sorted by its size
	Euler2D fluid;

	double gravity;
	//tolerance for level set
	double tol;

	//hash table to store particles
	HierarchicalHashTable hashedgrid;

	LevelSet combinedls;
	//Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal;

	Moving_RB() : particles(0), fluid(), gravity(0), combinedls() {} //, normal(0, 0) {}
	~Moving_RB() {};

	void add_sphere(const Domain2D&, const Coordinates&, double); //center and radius
	void add_particle(const Polygon&, const Domain2D&);
	void generate_hht(const Domain2D&); //generates the hierarchical hash table, also generating the list of springs (initialise to 0) for each particle

	int getNBodies(); //returns number of rigid body particles

	//-----------------------------------------------------
	//Construction of AABB (and?) sorting by size
	//-----------------------------------------------------
	void sort_by_size();

private:
	Moving_RB(const Moving_RB& rbsystem) : particles(rbsystem.particles), fluid(rbsystem.fluid), combinedls(rbsystem.combinedls) {}

};

#endif /* PARTICLES_H_ */