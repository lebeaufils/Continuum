#include "../headerFiles/Particles.h"

typedef std::unordered_multimap<std::pair<int, int>, int, HashPair> hmultimap;

//--------------------------------------------------------------
// Hashing
//--------------------------------------------------------------
//HashEntry
HashEntry::HashEntry(std::pair<int, int> key_input, int particle) : key(key_input) {
	labels.push_back(particle);
}

void HashEntry::add_label(int particle){
	labels.push_back(particle);
}

void HashEntry::remove_label(int particle){
	std::list<int>::iterator i;
	for (i=labels.begin(); i!=labels.end(); i++){
		if (*i == particle) break;
	}
	if (i != labels.end()){
		labels.erase(i);
	}
	//else - error check?
}

std::list<int> HashEntry::return_all_labels(int particle){
	std::list<int> label_list;
	std::list<int>::iterator i;
	for (i=labels.begin(); i!=labels.end(); i++){
		if (*i != particle){
			label_list.push_back(*i);
		}
	}
	return label_list;
}

//Hashtable
HashTable::HashTable(int res, double res_size) : resolution(res), resolution_size(res_size) {
	map = hmultimap();
	//--old--
	//buckets = make_unique<std::list<HashEntry>[]>(size); //creates an array of lists
}
/*
void HashTable::insert(std::pair< <int, int> > key, int particle){
	int index = hash_func(key); //hash the key 
	//check if the key already exists in the bucket
	std::list<HashEntry>::iterator i;
	for (i=buckets[index].begin(); i!=buckets[index].end(); i++){ //if list is empty, begin() == end()
		if (i->key == key) break; //increment the iterator until the key is found (or not)
	}
	//if the key is found in its indexed position
	if (i != buckets[index].end()){
		//append the particle label to the key's list
		i->add_label(particle);
	}
	else {
		//If the key dosent already exist, add the key to the hashed location
		buckets[index].pushback(HashEntry(key, particle)); //insert a new element in the linked list at the hashed index
	}
}

void HashTable::delete(std::pair< <int, int> > key, int particle){
	//find the hashed location of the given key
	int index = hash_func(key);
	//iterate through the linked list at the hashed index
	//to find the target key
	std::list<HashEntry>::iterator i;
	for (i=buckets[index].begin(); i!=buckets[index].end(); i++){
		if (i->key == key) break; //increment the iterator until the key is found (or not)
	}
	//if the key is found in its stored index, delete it from the linked list
	if (i != buckets[index].end()) {
		buckets[index].erase(i);
	}
	//error checking - delete should only be called knowing the particle exists in the grid
	//if it is not found, there is likely to be an error
	else {
		std::cout << "Grid (" << key.first << ", " << key.second << ") is empty.";	
	}
}

std::list<int> HashTable::search(std::pair< <int, int> > key, int particle){
	//find the hashed location of the given key
	int index = hash_func(key);
	//iterate through the linked list at the hashed index
	//to find the target key
	std::list<HashEntry>::iterator i;
	for (i=buckets[index].begin(); i!=buckets[index].end(); i++){
		if (i->key == key) break; //increment the iterator until the key is found (or not)
	}
	//if the key is found in its stored index
	if (i != buckets[index].end()) {
		i->return_all_labels(particle);
	}
	//error checking required
	else {
		std::cout << "Grid (" << key.first << ", " << key.second << ") is empty.";	
	}
}

int HashTable::hash_func(std::pair<int, int> key){
	boost::hash< std::pair<int, int> > hash; //to be replaced by perfect hash function
	return hash(key)%bucketsize;
}
*/

//Hierarchical hash table
HierarchicalHashTable::HierarchicalHashTable(int n) : tables(), close_tracker(), wall_tracker() { //where n is the number of particles in the system
	std::vector<int> particle_row(n, 0);
	std::vector<int> wall_row(2, 0); //2 represents horizontal or vertical boundary
	for (int i=0; i<n; i++){
		close_tracker.push_back(particle_row);
		wall_tracker.push_back(wall_row);
	}
}
//alternatively, if the default constructor was used
void HierarchicalHashTable::create_close_tracker(int n){
	std::vector<int> particle_row(n, 0);
	std::vector<int> wall_row(2, 0);
	for (int i=0; i<n; i++){
		close_tracker.push_back(particle_row);
		wall_tracker.push_back(wall_row);
	}
	//note: while the close tracker is a two dimensional array,
	//Since it is symmetric, it shall only be updated in the lower triangle format
}

void HierarchicalHashTable::add_table(double resolution_size){
	int resolution = static_cast<int>(tables.size());
	tables.push_back(HashTable(resolution, resolution_size));
	//adds an empty table to the hierarchy
}

//generates table based on particle sizes 
void HierarchicalHashTable::generate_tables(const Domain2D& domain, std::vector<Particle>& particlelist){
//takes a sorted list of particles
	std::vector<Particle>::iterator i = particlelist.begin();
	//generate a table for the largest particle, with resolution such that the particle occupies
	// between 0.5-1 of the grid cell: 0.5 <= sqrt(2*size*size)/sqrt(2*l*l) <= 1 assuming square grid
	//this reduces to. 0.5 <=  size/l <=1, where size is the length of the bounding box
	int label = 0; //index of particle in the particle list
	int resolution = 0; //resolution is zero for the largest particle
	double l = 1.5*i->size; //resolution length
	add_table(l);
	//Find the initial intersecting grid cells and store this information
	i->resolution = resolution;
	i->label = label;
	i->igrids = find_intersecting_cells(i->size, i->centre, resolution); 
	add_particle_all(domain, particlelist, *i, resolution);
	i++;
	label++;
	while (i!=particlelist.end()){
		//since the particles are decreasing in size, we only need to check for the lower bound
		if (i->size <= 0.5*l){ //if the particle occupies less than half of the grid cell,
			resolution++; //increase the resolution count
			l = 0.5*l; //halve the grid/bucket length
			add_table(l); //create a new layet on hash table
		}
		i->resolution = resolution;
		i->label = label;
		i->igrids = find_intersecting_cells(i->size, i->centre, resolution);
		add_particle_all(domain, particlelist, *i, resolution);
		i++;
		label++;
	}
	//This function adds a new layer to the Hierarchy with half the resolution, if the characteristic length of a particle becomes less than half the grid length
}

int HierarchicalHashTable::compute_resolution(double boxsize){
	//computes min resolution of particle
	int res = std::ceil(std::log(boxsize/tables[0].resolution_size)/std::log(2));
	return (res >= 0 ? res : 0);
}

std::vector<std::pair<int, int> > HierarchicalHashTable::find_intersecting_cells(double boxsize, const vector2& centre, int resolution){
	//calculate min(x,y) and max(x,y) of the bounding box
	vector2 minxy = centre - vector2(boxsize/2., boxsize/2.);
	vector2 maxxy = centre + vector2(boxsize/2., boxsize/2.);
	//Finds corners
	std::vector<std::pair<int, int> > grids;
	//corners.push_back(std::make_pair<int, int>(floor(minxy(0)/tables[resolution].resolution_size), floor(minxy(1)/tables[resolution].resolution_size))); //bottom left
	//corners.push_back(std::make_pair<int, int>(floor(maxxy(0)/tables[resolution].resolution_size), floor(maxxy(1)/tables[resolution].resolution_size))); //top right
	int iminus = floor(minxy(0)/tables[resolution].resolution_size);
	int iplus = floor(maxxy(0)/tables[resolution].resolution_size);
	int jminus = floor(minxy(1)/tables[resolution].resolution_size);
	int jplus = floor(maxxy(1)/tables[resolution].resolution_size);

	//grids to be added
	for (int i=iminus; i<=iplus; i++){
		for (int j=jminus; j<=jplus; j++){
			grids.push_back(std::make_pair(i, j)); //adds corner cells and those "between"
		}
	} 
	//The order of grid indices in <grids> are now bottom left, "sandwich cells", top right
	return grids;
}

void HierarchicalHashTable::close_to_wall(const Domain2D& domain, const Particle& gr, int resolution){
//To check for closeness to walls, consider the grids stored at each particle's respective resolutions
	//if the grid indices are at the edges, target particle is considered close to the walls and the counter is incremented
	int x_max = floor(domain.Lx/tables[resolution].resolution_size);
	int y_max = floor(domain.Ly/tables[resolution].resolution_size);

	std::vector<std::pair<int, int> >::const_iterator it;
	for (it=gr.igrids.begin(); it!=gr.igrids.end(); it++){
		//---- x axis
		if (it->first <= 0 || it->first >= x_max){
			wall_tracker[gr.label][0] = 1;
		}
		else { //reset to 0
			wall_tracker[gr.label][0] = 0;
		}
		//---- y axis
		if (it->second <= 0 || it->second >= y_max){
			wall_tracker[gr.label][1] = 1;
		}
		else { //reset
			wall_tracker[gr.label][1] = 0;
		}
	}
}

//add or remove particles
void HierarchicalHashTable::add_particle_all(const Domain2D& domain, const std::vector<Particle>& particlelist, const Particle& gr, int resolution){
//adds particle label to all cells for equal or lower resolution
//to be used in the table initialisation
	for (std::vector<std::pair<int, int> >::const_iterator i = gr.igrids.begin(); i!=gr.igrids.end(); i++){
		//check if other particles are already stored in the gridcell
		std::pair<hmultimap::iterator, hmultimap::iterator> itpair 
			= tables[resolution].map.equal_range(*i); //equal range returns the start and end iterators to pairs with the same key
		for (hmultimap::iterator it=itpair.first; it!=itpair.second; it++){
			//if the key (grid cell) where the particle is to be added already has other particles,
			//update the close tracker (since this resolution is res(gr))
			close_tracker[gr.label][it->second]++; //closetracker[a][b] where a is the row (particle being considered) and b is the col (particle being compared with)
			//close_tracker[it->second][gr.label]++; //unnecessary if we only consider the lower triangle
			//close_tracker has rows representing each particles' closeness tracker to every other particle
		}
		tables[resolution].map.insert({*i, gr.label}); //insert the label into relevant cells
	}
	//update wall tracker
	close_to_wall(domain, gr, resolution);

	//Need to check, before inserting, if there are other particles in the stored in the same grid.
	//if so, update the close_tracker, then insert
	//the same needs to be done for lower resolutions

	resolution--;
	while (resolution>0){
		//convert grid index to that of lower resolution
		//double res_ratio = tables[resolution].resolution_size/tables[resolution+1].resolution_size; //ratio of lengths of cells resolution-1/resolution
		double res_ratio = tables[resolution+1].resolution_size/tables[resolution].resolution_size;
		int iminus = floor(gr.igrids.front().first*res_ratio); 
		int iplus = floor(gr.igrids.back().first*res_ratio);
		int jminus = floor(gr.igrids.front().second*res_ratio);
		int jplus = floor(gr.igrids.back().second*res_ratio);

		//grids to be added
		for (int i=iminus; i<=iplus; i++){
			for (int j=jminus; j<=jplus; j++){
				//again, check for particles stored in the same grid key
				std::pair<hmultimap::iterator, hmultimap::iterator> itpair 
					= tables[resolution].map.equal_range({i, j});	
				for (hmultimap::iterator it=itpair.first; it!=itpair.second; it++){
					//if the key (grid cell) where the particle is to be added already has other particles,
					//check if the resolution is min(res(a), res(b)) where a and b are the two particles being tested
					if (resolution == fmin(particlelist[it->second].resolution, gr.resolution)){
						//if so, update the close tracker
						close_tracker[gr.label][it->second]++;
						//close_tracker[it->second][gr.label]++;
					}
					//close_tracker has rows representing each particles' closeness tracker to every other particle
				}
				
				tables[resolution].map.insert({{i, j}, gr.label});
			}
		} 
		resolution--;
	}

}

void HierarchicalHashTable::move_particle(const Domain2D& domain, const std::vector<Particle>& particlelist, Particle& gr){
//checks if the intersecting grid cells have changed
//assumes the particle's centre of mass has been updated
//New intersecting grids are computed using the updated centre and size of bounding box
	std::vector<std::pair<int, int> > newgrids = find_intersecting_cells(gr.size, gr.centre, gr.resolution);
	//note that the construction of the array of grid cells arre such that they are already sorted based on their i values, followed by j (i, j)
	//this follows the comparison rules for std::pair<int, int>
//we can thus find the intersection of the old and new list of grid cells by merging them
	std::vector<std::pair<int, int> >::const_iterator i = gr.igrids.begin(); //old intersecting grid
	std::vector<std::pair<int, int> >::const_iterator j = newgrids.begin(); //new intersecting grid

	//complexity of merging is O(m+n)
	std::vector<std::pair<int, int> > unchangedcells; //cells that remain unchanged after movemnet
	while(i!=gr.igrids.end() && j!=newgrids.end()){
		if (*i < *j) i++;
		else if (*i > *j) j++;
		else { //if i and j are the same element
			unchangedcells.push_back(*i);
			i++;
			j++;
		} //the intersection of both arrays is also a sorted array
	}
	//now that we have the intersection of both grid cells, we can find the cells where we need to add/ remove the particle label using set_difference
	//set difference returns the elements in v1 that are not found in v2 (non commutative)
	std::vector<std::pair<int, int> > remove_gr; //list of grids where gr is to be removed from the bucket
	std::vector<std::pair<int, int> > add_gr; //list of grids where gr is to be added to the bucket
	std::set_difference(gr.igrids.begin(), gr.igrids.end(), unchangedcells.begin(), unchangedcells.end(), std::inserter(remove_gr, remove_gr.begin()));
	std::set_difference(newgrids.begin(), newgrids.end(), unchangedcells.begin(), unchangedcells.end(), std::inserter(add_gr, add_gr.begin()));
	//remove or add label from the relevant buckets and update the close tracker
	int resolution = gr.resolution;
	remove_particle(particlelist, gr, resolution, remove_gr);
	add_particle(particlelist, gr, resolution, add_gr);

/////////////	//update the buckets of lower resolution, if applicable
	resolution--;
	while (resolution>0){ //check this
		//double res_ratio = tables[resolution].resolution_size/tables[resolution+1].resolution_size;
		double res_ratio = tables[resolution+1].resolution_size/tables[resolution].resolution_size;
		//reduce the removed grid keys to that of the lower resolution
		std::vector<std::pair<int, int> > lowres_remove;
		for (std::vector<std::pair<int, int> >::iterator it=remove_gr.begin(); it!=remove_gr.end(); it++){
			std::pair<int, int> index = {floor(it->first*res_ratio), floor(it->second*res_ratio)};
			//to prevent duplicate grid indices, append the new index only if the vector is empty, of it it is not the same as the last element
			if (lowres_remove.empty() || index != lowres_remove.back()){ 
				lowres_remove.push_back(index);
			}
		}
		//repeat for the keys to be added
		std::vector<std::pair<int, int> > lowres_add;
		for (std::vector<std::pair<int, int> >::iterator it=add_gr.begin(); it!=add_gr.end(); it++){
			std::pair<int, int> index = {floor(it->first*res_ratio), floor(it->second*res_ratio)};
			if (lowres_add.empty() || index != lowres_add.back()){ 
				lowres_add.push_back(index);
			}
		}
		//perform partical removal/addition with the new array of keys
		remove_particle(particlelist, gr, resolution, lowres_remove);
		add_particle(particlelist, gr, resolution, lowres_add);
		resolution--; //further lower the resolution
	}
	//end of while loop
	//update gr.igrids with the new intersecting grids
	gr.igrids = newgrids;
	close_to_wall(domain, gr, gr.resolution);
}

void HierarchicalHashTable::remove_particle(const std::vector<Particle>& particlelist, const Particle& gr, int resolution, const std::vector<std::pair<int, int> >& remove_gr){
	//removes particle label from keys stored in remove_gr, for the specified resolution
	for (std::vector<std::pair<int, int> >::const_iterator iterator=remove_gr.begin(); iterator!=remove_gr.end(); iterator++){
		//at each key, find the particle from all the elements stored in the grid and erase it
		std::pair<hmultimap::iterator, hmultimap::iterator> itpair = tables[resolution].map.equal_range(*iterator);
		for (hmultimap::iterator it=itpair.first; it!=itpair.second; it++){
			if (it!=itpair.second){
				if (it->second != gr.label){
					if (resolution == fmin(gr.resolution, particlelist[it->second].resolution)){
						//decrement the close tracker
						close_tracker[gr.label][it->second]--;
						//close_tracker[it->second][gr.label]--;
					}
					//if a close tracker is 0 -- there are no longer any buckets that contain both labels
				}
				//else tables[resolution].map.erase(it);
			}
			else {
				std::cout << "remove particle failed, Grid is empty" << std::endl;
			}
		}
	}

	//this removes the label from the actual resolution. what of tables of lower resolution?
	//can still call this function, simply lower the resolution and find the corresponding grid indices before calling.
	//this can be done by diving the grid index by the ratios of resolution length, flooring the result. Add to a new vector any pairs that are not duplicates
}

void HierarchicalHashTable::add_particle(const std::vector<Particle>& particlelist, const Particle& gr, int resolution, const std::vector<std::pair<int, int> >& add_gr){
	//adds particle labels to keys stored in addgr
	for (std::vector<std::pair<int, int> >::const_iterator iterator=add_gr.begin(); iterator!=add_gr.end(); iterator++){
		std::pair<hmultimap::iterator, hmultimap::iterator> itpair = tables[resolution].map.equal_range(*iterator);
		for (hmultimap::iterator it=itpair.first; it!=itpair.second; it++){
			if (it!=itpair.second){ //checks if the key is not empty
				if (resolution == fmin(gr.resolution, particlelist[it->second].resolution)){
					close_tracker[gr.label][it->second]++;
					//close_tracker[it->second][gr.label]++;
				}
			}
		}
		//for each key, insert the label into the hash map
		tables[resolution].map.insert({*iterator, gr.label});
	}	
}

//--------------------------------------------------------------
//Particles
//--------------------------------------------------------------
//circle constructor
Particle::Particle(const Domain2D& domain, const Coordinates& center, double r) : mass(1.0), density(1.0), damping_coefficient(0.0), 
miu(0.3), k_n(50), k_s(50), k_c(50), label(0), size(0), resolution(0), igrids(0), ls(), centroid(0, 0), ref_nodes(0), 
dynamicls(), centre(0, 0), nodes(0), vc(0, 0), w(0), s(0, 0), theta(0), force(0, 0), springs(), wall_springs(), overlap(), wall_overlap(), torque(0) {
	LevelSetMethods::initialise_circle(ls, domain, center.x, center.y, r);
	dynamicls = ls;
	centroid = center_of_mass(ls, *this, domain);
	centre = centroid;

	//seeding nodes for the circle
	//double circumference = 2*r*atan(1.0)*4;
	//double spacing = 2*r/10.;
	//vector2 surface_p(center.x + r, center.y);

	vector2 surface_p(center.x + r, center.y);
	//find the first point on the zeroth levelset contour
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (ls.phi(i+domain.buffer, j+domain.buffer) == 0){
				surface_p(0) = domain.X(i, j).x;
				surface_p(1) = domain.X(i, j).y;
				break;
			}
		}
	}
	ref_nodes.push_back(surface_p);
	nodes.push_back(surface_p);
	//find the frequency of rotation
	int nodesize = 32; //circumference/spacing = 10pi
	double t = 0;
	for (int dt=1; dt<nodesize; dt++){
		vector2 v = Rotor2::rotate_about(surface_p, vector2(center.x, center.y), 2*3.14159/nodesize, t+dt); 
		ref_nodes.push_back(v);
		nodes.push_back(v);
	}

	//calculate the size of the bounding box
	size = AABB_extent(*this);		

	//initialise the horizontal and vertical wall springs to 0
	std::unordered_map<int, vector2> nodesprings;
	std::unordered_map<int, double> nodedist;
	wall_springs.push_back(nodesprings);
	wall_springs.push_back(nodesprings);
	wall_overlap.push_back(nodedist);
	wall_overlap.push_back(nodedist);
}
//polygon constructor
Particle::Particle(const Polygon& poly, const Domain2D& domain) : mass(1.0), density(1.0), damping_coefficient(0.2), 
miu(0.3), k_n(20), k_s(20), k_c(50), label(0), size(0), resolution(0), igrids(0), ls(), centroid(0, 0), ref_nodes(0), 
dynamicls(), centre(0, 0), nodes(0), vc(0, 0), w(0), s(0, 0), theta(0), force(0, 0), springs(), wall_springs(), overlap(), wall_overlap(), torque(0) {
	LevelSetMethods::initialise(ls, domain, poly);
	dynamicls = ls;
	centroid = center_of_mass(ls, *this, domain);
	centre = centroid;

	//double circumference = static_cast<int>(poly.surfacepoints.size())*domain.dx;
	//double diameter = circumference/3.2;
	//int spacing = floor(diameter/(10*domain.dx));
	int spacing = floor(static_cast<int>(poly.surfacepoints.size())/32);
	for (int a=0; a < static_cast<int>(poly.surfacepoints.size()); a+=spacing) {
		Coordinates tmp = domain.X(poly.surfacepoints[a].i, poly.surfacepoints[a].j);
		ref_nodes.push_back(vector2(tmp.x, tmp.y));
		nodes.push_back(vector2(tmp.x, tmp.y));
	}

	//calculate the size of the bounding box
	size = 1.5*AABB_extent(*this);	

	//initialise the horizontal and vertical wall springs to 0
	std::unordered_map<int, vector2> nodesprings;
	wall_springs.push_back(nodesprings);
	wall_springs.push_back(nodesprings);
	wall_overlap.push_back(nodedist);
	wall_overlap.push_back(nodedist);
}
//copy constructor
Particle::Particle(const Particle& gr) : density(gr.density), damping_coefficient(gr.damping_coefficient), 
miu(gr.miu), k_n(gr.k_n), k_s(gr.k_s), k_c(gr.k_c), label(gr.label), size(gr.size), resolution(gr.resolution), igrids(gr.igrids), ls(gr.ls), centroid(gr.centroid), ref_nodes(gr.ref_nodes), 
dynamicls(gr.dynamicls), centre(gr.centre), nodes(gr.nodes), vc(gr.vc), w(gr.w), s(gr.s), theta(gr.theta), force(gr.force), springs(gr.springs), 
wall_springs(gr.wall_springs), overlap(gr.overlap), wall_overlap(gr.wall_overlap), torque(gr.torque) {
	 for (int a=0; a<static_cast<int>(gr.nodes.size()); a++){
	 	nodes.push_back(gr.nodes[a]);
	 	ref_nodes.push_back(gr.ref_nodes[a]);
	 }
}

void Particle::set_mass(double m){
	mass = m;
}

void Particle::set_density(double d){
	density = d;
}

void Particle::set_velocity(const vector2& trans_velocity, double angular_velocity){
	vc = trans_velocity;
	w = angular_velocity;
}

double Particle::compute_mass(const Particle& gr, const Domain2D& domain){

	double m = 0;

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			m += LevelSetMethods::smoothed_heaviside(gr.ls, domain, i+domain.buffer, j+domain.buffer);
			//sums the number of cells within the rigid body, including the
			//transition zone
		}
	}
	m = gr.density*domain.dx*domain.dy*m; //rho * g^2 where g is the grid spacing

	return m;
}

vector2 Particle::center_of_mass(const LevelSet& ls, const Particle& gr, const Domain2D& domain){
	//For a mass with density rho(r) within the solid, the integral of the weighted (density) position
	//coordinates relative to its center of mass is 0
	double c_x = 0;
	double c_y = 0;

	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			c_x += LevelSetMethods::smoothed_heaviside(ls, domain, i+domain.buffer, j+domain.buffer)*(domain.X(i, j).x);
			c_y += LevelSetMethods::smoothed_heaviside(ls, domain, i+domain.buffer, j+domain.buffer)*(domain.X(i, j).y);
		}
	}
	double m = compute_mass(gr, domain);
	c_x = (gr.density*domain.dx*domain.dy/m)*c_x;
	c_y = (gr.density*domain.dx*domain.dy/m)*c_y;
	vector2 c(c_x, c_y);
	return c;
}

double Particle::moment_of_inertia(const LevelSet& ls, const Particle& gr, const Domain2D& domain){
	//for rotation confined to a plane, the mass moment of inertia reduces to a scalar value.
	vector2 c = center_of_mass(ls, gr, domain);
	//in 2d, moment of inertia is a 1x1 tensor given by J = sum(m_i * (x_i^2 + y_i^2))
	double inertia = 0;
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			inertia += LevelSetMethods::smoothed_heaviside(ls, domain, i+domain.buffer, j+domain.buffer)*(pow((i*domain.dx - c(0)),2) + pow((j*domain.dy - c(1)),2));
		}
	}
	inertia = -gr.density*domain.dx*domain.dy*inertia;
	return inertia;
}

vector2 Particle::velocity(const Coordinates& p, const Particle& gr){
	vector2 v(0, 0);
	vector2 r(p.x-gr.centroid(0), p.y-gr.centroid(1));
	//v = vc + w cross r
	v(0) = gr.vc(0) - r(1)*gr.w;
	v(1) = gr.vc(1) + r(0)*gr.w;
	return v;
}

LevelSet Particle::merge(const std::vector<Particle>& particles, const Domain2D& domain){
	LevelSet ls;
	//initialise the levelset with a large positive number
	ls.phi = Eigen::MatrixXd::Constant(domain.Nx+2*domain.buffer, domain.Ny+2*domain.buffer, 1e6);

	for (int i=0; i<domain.Nx+2*domain.buffer; i++){
		for (int j=0; j<domain.Ny+2*domain.buffer; j++){
			double min_phi = 1e6;
			for (int a=0; a<static_cast<int>(particles.size()); a++){
				if (particles[a].dynamicls.phi(i, j) < min_phi) min_phi = particles[a].dynamicls.phi(i, j);
			}
			ls.phi(i, j) = min_phi;
		}
	}
	//LevelSetMethods::boundary_conditions(ls, domain);

	return ls;
}

vector2 Particle::normal_sum(const LevelSet& ls, const Domain2D& domain){
	vector2 n(0,0);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			vector2 n_i = LevelSetMethods::normal(ls, domain, i+domain.buffer, j+domain.buffer);
			double delta = LevelSetMethods::smoothed_delta(ls, domain, i+domain.buffer, j+domain.buffer);
			//std::cout << delta << '\t';
			n += n_i*delta*domain.dx*domain.dy;
		}
		//std::cout << std::endl;
	}
	return n;
}

double Particle::length_sum(const LevelSet& ls, const Domain2D& domain){
	double s=0;
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			double delta = LevelSetMethods::smoothed_delta(ls, domain, i+domain.buffer, j+domain.buffer);
			//std::cout << delta << '\t';
			s += domain.dx*domain.dy*delta;
		}
		//std::cout << std::endl;
	}
	return s;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
LevelSet Particle::motion(const Domain2D& domain, const vector2& s, double theta){
	//Temporary storage for levelset values at time t
	//The levelset can be updated using its initial values and current position
	LevelSet ls_t; //ls at time t
	ls_t.phi = Eigen::MatrixXd::Constant(domain.Nx+2*domain.buffer, domain.Ny+2*domain.buffer, 1e6);//matrix::Zero(domain.Nx+2*domain.buffer, domain.Ny+2*domain.buffer);

	vector2 newcentre = centroid + s;
	//calculate levelset within bounding box
	int min_i = floor((newcentre(0)-size/1.8)/domain.dx);
	int max_i = ceil((newcentre(0)+size/1.8)/domain.dx);
	int min_j = floor((newcentre(1)-size/1.8)/domain.dy);
	int max_j = ceil((newcentre(1)+size/1.8)/domain.dy);


	for (int i=min_i; i<=max_i; i++){
		for (int j=min_j; j<=max_j; j++){
			Coordinates originalpos(i*domain.dx, j*domain.dy);
			originalpos = LevelSetMethods::translation_reverse(originalpos, s); //translate the point
			originalpos = LevelSetMethods::rotation_reverse(originalpos, centroid, theta); //rotate the point
			ls_t.phi(i+domain.buffer, j+domain.buffer) = LevelSetMethods::interpolation_value(ls, domain, originalpos); //interpolate value from reference ls
		}	
	}
	
	//move the nodes from ref, store a seperate node list
	for (int a=0; a<static_cast<int>(nodes.size()); a++){
		//std::cout << nodes[a].transpose() << std::endl;
		Coordinates node_a(ref_nodes[a]); //always start from reference nodes
		node_a = LevelSetMethods::rotation(node_a, centroid, theta);
		node_a = LevelSetMethods::translation(node_a, s); 
		nodes[a] = vector2(node_a.x, node_a.y);
	}
	
////////
	LevelSetMethods::fast_sweep(ls_t, domain);
	//LevelSetMethods::boundary_conditions(ls_t, domain);
	return ls_t;
}

vector2 Particle::cross(double w, const vector2& pos){
 	return vector2(-w*pos(1), w*pos(0));
}

double Particle::cross(const vector2& m, const vector2& n){
	return m(0)*n(1) - m(1)*n(0);
}

double Particle::AABB_extent(const Particle& gr){
	double min_x=1e6;
	double max_x=0;
	double min_y=1e6;
	double max_y=0;
	
	/*for (int i=0; i<ls.phi.rows(); i++){
		for (int j=0; j<ls.phi.cols(); j++){
			if (ls.phi(i, j) < size) size = ls.phi(i, j);
		}
	}*/
	//using the nodes to find bounding box
	for (int i=0; i<static_cast<int>(gr.ref_nodes.size()); i++){
		if (gr.ref_nodes[i](0) > max_x) max_x = gr.ref_nodes[i](0);
		if (gr.ref_nodes[i](0) < min_x) min_x = gr.ref_nodes[i](0);
		if (gr.ref_nodes[i](1) > max_y) max_y = gr.ref_nodes[i](1);
		if (gr.ref_nodes[i](1) < min_y) min_y = gr.ref_nodes[i](1);
	}

	//currently, size measure of the furthest point of the particle from its centroid.
	//double size to get the bounding box.
	double size = fmax(max_x - min_x, max_y - min_y);
	return size; 
}

///////
void Moving_RB::add_particle(const Polygon& poly, const Domain2D& domain){
	//Particle new_particle(poly, domain);
	particles.push_back(Particle(poly, domain));
}

void Moving_RB::add_sphere(const Domain2D& domain, const Coordinates& center , double r){
	//Particle new_particle(domain, center, r);
	particles.push_back(Particle(domain, center, r));
}

void Moving_RB::sort_by_size(){
	//sorts the particles based on tHier bounding box(square) size
	std::sort(particles.begin(), particles.end(), [](const Particle& a, const Particle& b) {
        return a.size > b.size;   
    });
}

//After all particles have been added
void Moving_RB::generate_hht(const Domain2D& domain){ //generate hierarchical hash table
	int n = static_cast<int>(particles.size());
	hashedgrid.create_close_tracker(n);
	hashedgrid.generate_tables(domain, particles);

	//particle spring
	std::vector<Particle>::iterator it;
	for (it = particles.begin(); it != particles.end(); it++){
		for (int i=0; i<n; i++){
			std::unordered_map<int, vector2> nodesprings;
			std::unordered_map<int, double> nodedist;
			it->springs.push_back(nodesprings);
			it->overlap.push_back(nodedist);
		}
	}
}
//during contact collision, for each particle loop through the close counter and only perform the collision detection on particles > 0 in the close counter
//if there is contact, calculate the tangential spring and store it in the array under the respective label.


