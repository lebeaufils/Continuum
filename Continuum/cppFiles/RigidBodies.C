#include "../headerFiles/RigidBodies.h"

void RigidBodies::boundary_conditions(vecarray &U, const Domain2D& domain){

	if (U.cols() == 0 || U.rows() == 0){
		throw "Array is empty.";
	}

	//assigning ghost values in the x-direction 
	for (int j=0; j<domain.Ny; j++){
		U(1, j+2) = U(2, j+2);
		U(0, j+2) = U(1, j+2);
		U(domain.Nx+2, j+2) = U(domain.Nx+1, j+2);
		U(domain.Nx+3, j+2) = U(domain.Nx+2, j+2);
	} 
	//assigning ghost values in the y-direction
	for (int i=0; i<domain.Nx; i++){
		U(i+2, 1) = U(i+2, 2);
		U(i+2, 0) = U(i+2, 1);
		U(i+2, domain.Ny+2) = U(i+2, domain.Ny+1);
		U(i+2, domain.Ny+3) = U(i+2, domain.Nx+2);
	} 
}

void RigidBodies::reflected_state(Stationary_RB &var, const vecarray &primV, int i, int j, const vector2 &n_i){
	vector2 velocity(primV(i,j)(1), primV(i,j)(2));
	vector2 v_ghost = velocity - 2*(n_i.dot(velocity))*(n_i);

	vector4 reflected(primV(i,j)(0), v_ghost(0), v_ghost(1), primV(i,j)(3));
	vector4 reflected_cons = var.fluid.state_function->conservedVar2Dx(reflected);

	var.fluid.U(i, j) = reflected_cons;
}

//could compress this using functors, temporarily creating duplicate function
void RigidBodies::fast_sweep(const LevelSet &ls, Stationary_RB &rb, const Domain2D &domain, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> &normal){

	//Multi dimensional extrapolation using fast sweeping
	auto heaviside = [ls](int i, int j){
		if (ls.phi(i, j) < 0) return 1;
		else return 0;
	};

	auto ghost_fluid = [domain](vecarray &primitive, int i, int j, vector2 n_i){
		//info is swept from outside into the rigid body
		//constant extrapolation of density (d), pressure (P) and velocity (v_i)
		//using the first order upwind finite difference approximation
		double dij=0, Pij=0, uij=0, vij=0;
		double dx, Px=0, ux=0, vx=0;
		double dy, Py=0, uy=0, vy=0;
		//consider the x derivatives
		if (n_i(0) > 0){
			//Information travels from the fluid to the rigid body (+ve to -ve levelset)
			//information is travelling to the left if the the x component of the normal is > 0
			//as this indicates that the rigid body is left of the fluid
			//right side is upwind
			dx = primitive(i+1,j)(0);
			ux = primitive(i+1,j)(1);
			vx = primitive(i+1,j)(2);
			Px = primitive(i+1,j)(3);
		}
		else {
			//information is travelling to the right
			//left side is upwind
			dx = primitive(i-1,j)(0);
			ux = primitive(i-1,j)(1);
			vx = primitive(i-1,j)(2);
			Px = primitive(i-1,j)(3);
		}

		if (n_i(1) > 0){
			dy = primitive(i,j+1)(0);
			uy = primitive(i,j+1)(1);
			vy = primitive(i,j+1)(2);
			Py = primitive(i,j+1)(3);
		}

		else {
			dy = primitive(i,j-1)(0);
			uy = primitive(i,j-1)(1);
			vy = primitive(i,j-1)(2);
			Py = primitive(i,j-1)(3);
		}
		//for constant extrapolation f(i, j) = 0
		double nx = n_i(0), ny = n_i(1);
		if (dx == 0) nx = 0;
		if (dy == 0) ny = 0;
		//using the update formula for q(i, j)
		double nxny = 0;
		if ((nx!=0) || (ny!=0)){
			if((nx>0) != (ny>0)){
				nxny = nx/domain.dx - ny/domain.dy;
				dij = (nx*(dx)/domain.dx - ny*(dy)/domain.dy)/nxny;
				Pij = (nx*(Px)/domain.dx - ny*(Py)/domain.dy)/nxny;
				uij = (nx*(ux)/domain.dx - ny*(uy)/domain.dy)/nxny;
				vij = (nx*(vx)/domain.dx - ny*(vy)/domain.dy)/nxny;
			}
			else {
				nxny = nx/domain.dx + ny/domain.dy;
				dij = (nx*(dx)/domain.dx + ny*(dy)/domain.dy)/nxny;//(n_i(0)*(dx)/domain.dx + n_i(1)*(dy)/domain.dy)/nxny;
				Pij = (nx*(Px)/domain.dx + ny*(Py)/domain.dy)/nxny;
				uij = (nx*(ux)/domain.dx + ny*(uy)/domain.dy)/nxny;
				vij = (nx*(vx)/domain.dx + ny*(vy)/domain.dy)/nxny;
			}
		}

		//std::cout << dij << '\t' << i-2 << ", " << j-2 << '\t' << '\t' << dx << '\t' << dy << '\t' << '\t' << nxny << '\t' << nx << '\t' << ny << std::endl; //<< Pij << '\t' << uij << '\t' << vij << std::endl;

		primitive(i, j) = vector4(dij, uij, vij, Pij);
	};

	//Convert the conserved variables into primitive form
	//vecarray primitive(domain.Nx+4, domain.Ny+4);
	//for (int i=0; i<domain.Nx+4; i++){
	//	for (int j=0; j<domain.Ny+4; j++){
	//		primitive(i, j) = rb.fluid.state_function->primitiveVar(rb.fluid.U(i, j));
	//	}
	//}
	//Convert the conserved variables into primitive form
	vecarray primitive(domain.Nx+4, domain.Ny+4);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			primitive(i+2, j+2) = rb.fluid.state_function->primitiveVar(rb.fluid.U(i+2, j+2));
		}
	}

	boundary_conditions(primitive, domain);

	//Gauss Seidel iteration for alternate directions.
	//A total of four sweeps is performed over the entire computational domain
	//1) i = 1:I, j = 1:J
	//2) i = I:1, j = 1:J
	//3) i = I:1, J = J:1
	//4) i = 1:I, j = J:1

	// a narrowband of abs(phi) < cellwidth is imposed
	auto gauss_seidel = [heaviside, ghost_fluid, domain, ls, normal](vecarray &primitive){
		//1) i = 1:I, j = 1:J
		for (int i=0; i<domain.Nx; i++){
			for (int j=0; j<domain.Ny; j++){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//2) i = I:1, j = 1:J
		for (int i=domain.Nx-1; i>=0; i--){
			for (int j=0; j<domain.Ny; j++){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//3) i = I:1, J = J:1
		for (int i=domain.Nx-1; i>=0; i--){
			for (int j=domain.Ny-1; j>=0; j--){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//1) i = 1:I, j = 1:J
		for (int i=0; i<domain.Nx; i++){
			for (int j=domain.Ny-1; j>=0; j--){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}
	};

	gauss_seidel(primitive);
	//gauss_seidel(primitive);

//////Implement narrow band here
	//reflect the ghostfluid states
	for (int i=0; i<domain.Nx; i++){
		for(int j=0; j<domain.Ny; j++){
			if (ls.phi(i+1, j+1) < 0) //&& ls.phi(i, j) > -(domain.dx + domain.dy)) //&& ls.phi(i+1, j+1) >= -(Test.domain.dx))
				reflected_state(rb, primitive, i+2, j+2, normal(i, j));
			}
	}
}

void RigidBodies::initial_conditions(rigidTests &Test){
	//Setting the initial fluid conditions
	Test.var.fluid.U.resize(Test.domain.Nx+4, Test.domain.Ny+4);
	Test.var.fluid.F.resize(Test.domain.Nx+2, Test.domain.Ny+2);
	Test.var.fluid.G.resize(Test.domain.Nx+2, Test.domain.Ny+2);

	//When solving for the fluid values, we can consider all the levelsets at once by taking their union
	Test.var.combinedls = LevelSetMethods::merge(Test.var.levelsets, Test.domain);
//////////////needs optimising/////////////
	for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			//for (int a=0; a<static_cast<int>(Test.var.levelsets.size()); a++){ 
				if (Test.var.combinedls.phi(i+1, j+1) >= 0){
					if (Test.interface(i, j) == false){
						Test.var.fluid.U(i+2, j+2) = Test.var.fluid.state_function->conservedVar2Dx(Test.initialL);
					}
					else {
						Test.var.fluid.U(i+2, j+2) = Test.var.fluid.state_function->conservedVar2Dx(Test.initialR);
					}
				}

				else {
					Test.var.fluid.U(i+2, j+2).setZero();
				}
			//}
			Test.var.fluid.F(i+1, j+1).setZero();
			Test.var.fluid.G(i+1, j+1).setZero();
		}
	}

	MUSCL::boundary_conditions(Test.var.fluid, Test.domain);

	//compute the normal vector using the levelset function
	Test.var.normal.resize(Test.domain.Nx, Test.domain.Ny);
	for (int i=0; i<Test.domain.Nx; i++){
		for (int j=0; j<Test.domain.Ny; j++){
			Test.var.normal(i, j) = LevelSetMethods::normal(Test.var.combinedls, Test.domain, i+1, j+1);
		}
	}
	//populate the rigid body with ghost values
	//note that these values need to be reflected before use
	//values are currently reflected in the fast sweep function
	fast_sweep(Test.var.combinedls, Test.var, Test.domain, Test.var.normal);
}

void RigidBodies::solver(Stationary_RB &var, Domain2D &domain, double CFL){
	//This is identical to the standard MUSCL solver
	//with the exception of ignoring points with levelset < 0
	//and potentially dropping to a first order scheme when the rigid body
	//involves a single grid point (MUSCL requires 1 real and 2 ghost points)
	
	//Storage for (piecewise linear) interpolated cell values 
	matrix ULix(domain.Nx+4, 4);
	matrix URix(domain.Nx+4, 4);

	matrix ULiy(domain.Ny+4, 4);
	matrix URiy(domain.Ny+4, 4);

	//Temporary storage for x and y conserved variables 
	matrix Uxn(domain.Nx+4, 4);
	matrix Uyn(domain.Ny+4, 4);


	slopeLimiter a = VanLeer;//getLimiter();

	double t = 0.0;
	int count = 0;
	do{
		//set timestep, taking maximum wavespeed in x or y directions
		//looping through the entire domain, find the maximum possible wavespeed in both directions
		double Smax=0;
		//double Smax_y=0;
		for (int i=2; i<domain.Nx+2; i++){
			for (int j=2; j<domain.Ny+2; j++){
				vector4 Utmp = var.fluid.U(i, j);
				double a_ij = var.fluid.state_function->soundspeed(Utmp);
				double d_ij_x = Utmp(1)/Utmp(0); //velocity in x direction
				double d_ij_y = Utmp(3)/Utmp(0); //velocity in y direction
				double vn = sqrt(pow(d_ij_x, 2) + pow(d_ij_y, 2));
				double lambda = fmax(fmax(abs(d_ij_x) + a_ij, abs(d_ij_y) + a_ij), vn + a_ij);
				if (lambda > Smax) Smax = lambda;
			}
		}
		if (count <= 5){
				//using a CFL number of 0,.2 for the first 5 timesteps
				domain.dt = 0.2*fmin((domain.dx/Smax), (domain.dy/Smax));
			}
		else {
				domain.dt = CFL*fmin((domain.dx/Smax), (domain.dy/Smax));
		}
		if (t + domain.dt > domain.tstop) domain.dt = domain.tstop - t;

		t += domain.dt;
		count += 1;

		//sweeping in the x-direction for each y row
		for (int j=0; j<domain.Ny; j++){
			//Within each row, store the values for each grid point in x with a single row list
			for(int i=0; i<domain.Nx+4; i++){
				Uxn.row(i) = var.fluid.U(i, j+2);
				//the sweep in x is only performed within the real y grid
			} 
			//calculate and update to new interpolated values U_i
			MUSCL::data_reconstruction(Uxn, a, ULix, URix, domain.Nx);

			//sweeping in the x-direction for each y row
			for (int i=1; i<domain.Nx+2; i++){
				if (var.combinedls.phi(i, j+1) >= 0 || var.combinedls.phi(i+1, j+1) >= 0 || var.combinedls.phi(i-1, j+1) >= 0){
					vector4 Fx(0, 0, 0, 0);
					//edge case where there is only one ghost cell
					if (var.combinedls.phi(i+1, j+1) < 0 && var.combinedls.phi(i+2, j+1) >= 0){
						MUSCL::compute_fluxes(var.fluid, Uxn, Fx, i);
						//std::cout << i-1 << '\t' << j << '\t' << var.combinedls.phi(i+1, j+1) << std::endl;
					}
					//regular case
					else {
						MUSCL::compute_fluxes(var.fluid, domain, Uxn, Fx, i, ULix, URix, domain.dx);
						//std::cout << Fx.transpose() << '\t' << '(' << i << ", " << j+1 << ')' << std::endl;
						//if (i==13) std::cout << Uxn.row(i) << '\t' << '\t' << ULix.row(i) << '\t' << '\t'  << URix.row(i) << std::endl;
					}
					var.fluid.F(i, j+1) = Fx; //storing the computed flux
				}
			}
		}


		//sweeping in the y-direction for each x row
		for (int i=0; i<domain.Nx; i++){
			for (int j=0; j<domain.Ny+4; j++){
				//Since the normal velocity is now v, its position is swapped with u
				Uyn.row(j) = var.fluid.swap_xy(var.fluid.U(i+2, j));
			}
			MUSCL::data_reconstruction(Uyn, a, ULiy, URiy, domain.Ny);

			for (int j=1; j<domain.Ny+2; j++){
				if (var.combinedls.phi(i+1, j) >= 0 || var.combinedls.phi(i+1, j+1) >= 0 || var.combinedls.phi(i+1, j-1) >=0){
					vector4 Fy(0, 0, 0, 0);
					if (var.combinedls.phi(i+1, j+1) < 0 && var.combinedls.phi(i+1, j+2) >= 0){
						MUSCL::compute_fluxes(var.fluid, Uyn, Fy, j);
						//std::cout << i << '\t' << j-1 << '\t' << var.combinedls.phi(i+1, j+1) << std::endl;
					}	
					else {
						MUSCL::compute_fluxes(var.fluid, domain, Uyn, Fy, j, ULiy, URiy, domain.dy);
					}
					var.fluid.G(i+1, j) = Fy; //storing the computed flux.
				}		
			}
		}

		//updating U with Strand splitting -- X(0.5dt)Y(dt)X(0.5dt)
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				if (var.combinedls.phi(i-1, j+1) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i, j+2), var.fluid.F(i, j+1), var.fluid.F(i-1, j+1), domain.dt/2, domain.dx);
					//std::cout << var.fluid.F(i, j+1).transpose()  << '\t'  << '\t' << var.fluid.F(i-1, j+1).transpose() << std::endl;
				}
			}
		}
		for (int i=0; i<domain.Nx; i++){
			for (int j=2; j<domain.Ny+2; j++){
				if (var.combinedls.phi(i+1, j-1) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i+2, j), var.fluid.swap_xy(var.fluid.G(i+1, j)), var.fluid.swap_xy(var.fluid.G(i+1, j-1)), domain.dt, domain.dy);
				}
			}
		}		
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				if (var.combinedls.phi(i-1, j+1) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i, j+2), var.fluid.F(i, j+1), var.fluid.F(i-1, j+1), domain.dt/2, domain.dx);
				}
			}
		}
		MUSCL::boundary_conditions(var.fluid, domain);

		//Extrapolate new ghost values
		fast_sweep(var.combinedls, var, domain, var.normal);

		if (count%50==0) std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << Smax_x << '\t' << Smax_y << std::endl;
		//if (count == 1) t = domain.tstop;
	}while (t < domain.tstop);
	std::cout << "count = "  << count << std::endl;
}

void RigidBodies::output(const Stationary_RB &var, const Domain2D &domain){

	std::ofstream outfile;
	outfile.open("dataeuler.txt");

	for (int i=2; i<domain.Nx+2; i++){
		for (int j=2; j<domain.Ny+2; j++){
			if (var.levelsets[0].phi(i-1, j-1) >= 0){
				vector4 Ux = var.fluid.U(i, j);
				double u = sqrt(pow(Ux(1)/Ux(0), 2) + pow(Ux(3)/Ux(0), 2));
				double P = var.fluid.state_function->Pressure(Ux);
				double e = var.fluid.state_function->internalE(Ux);

				//central difference to calculate partial derivatives in x, y for density
				double grad_density_x = (var.fluid.U(i+1, j)(0) - var.fluid.U(i-1, j)(0))/(2*domain.dx);
				double grad_density_y = (var.fluid.U(i, j+1)(0) - var.fluid.U(i, j-1)(0))/(2*domain.dy);
				//calculating the numerical schlieren
				double schlieren = exp((-20*sqrt(pow(grad_density_x, 2) + pow(grad_density_y, 2)))/(1000*Ux(0)));

				outfile << domain.dx*(i-2) << '\t' << domain.dy*(j-2) << '\t' << Ux(0) << '\t' << u
				<< '\t' << P << '\t' << e << '\t' << schlieren << '\t' << '\t' << '\t' << var.levelsets[0].phi(i-1, j-1) << std::endl;
			}
			/*else {
				outfile << domain.dx*(i-2) << '\t' << domain.dy*(j-2) << '\t' << 0 << '\t' << 0
				<< '\t' << 0 << '\t' << 0 << '\t' << 0 << std::endl;
			}*/
		}
		outfile << std::endl;
	}

	outfile.close();
	std::cout << "done: Rigid Body" << std::endl;
}

void RigidBodies::rigid_body_solver(rigidTests &Test, double CFL){
	initial_conditions(Test);
	solver(Test.var, Test.domain, CFL);
	output(Test.var, Test.domain);
}


//------------------------------------------------------------------------
// Moving Level Set and Forces
//------------------------------------------------------------------------
void RigidBodies::reflected_state(Moving_RB &var, const vecarray &primV, int i, int j, const vector2 &n_i, const vector2& vb){
	//Contribution from the moving body is considered
	vector2 velocity(primV(i,j)(1), primV(i,j)(2));
	vector2 v_ghost = velocity - 2*(n_i.dot(velocity))*(n_i) + 2*(n_i.dot(vb))*(n_i);

	vector4 reflected(primV(i,j)(0), v_ghost(0), v_ghost(1), primV(i,j)(3));
	vector4 reflected_cons = var.fluid.state_function->conservedVar2Dx(reflected);

	var.fluid.U(i, j) = reflected_cons;
}

void RigidBodies::fast_sweep(const LevelSet& ls, const Particle& gr, Moving_RB &rb, const Domain2D &domain, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> &normal){

	//Multi dimensional extrapolation using fast sweeping
	auto heaviside = [ls, gr](int i, int j){
		if (ls.phi(i, j) < 0) return 1;
		else return 0;
	};

	auto ghost_fluid = [domain](vecarray &primitive, int i, int j, vector2 n_i){
		//info is swept from outside into the rigid body
		//constant extrapolation of density (d), pressure (P) and velocity (v_i)
		//using the first order upwind finite difference approximation
		double dij=0, Pij=0, uij=0, vij=0;
		double dx, Px=0, ux=0, vx=0;
		double dy, Py=0, uy=0, vy=0;
		//consider the x derivatives
		if (n_i(0) > 0){
			//Information travels from the fluid to the rigid body (+ve to -ve levelset)
			//information is travelling to the left if the the x component of the normal is > 0
			//as this indicates that the rigid body is left of the fluid
			//right side is upwind
			dx = primitive(i+1,j)(0);
			ux = primitive(i+1,j)(1);
			vx = primitive(i+1,j)(2);
			Px = primitive(i+1,j)(3);
		}
		else {
			//information is travelling to the right
			//left side is upwind
			dx = primitive(i-1,j)(0);
			ux = primitive(i-1,j)(1);
			vx = primitive(i-1,j)(2);
			Px = primitive(i-1,j)(3);
		}

		if (n_i(1) > 0){
			dy = primitive(i,j+1)(0);
			uy = primitive(i,j+1)(1);
			vy = primitive(i,j+1)(2);
			Py = primitive(i,j+1)(3);
		}

		else {
			dy = primitive(i,j-1)(0);
			uy = primitive(i,j-1)(1);
			vy = primitive(i,j-1)(2);
			Py = primitive(i,j-1)(3);
		}
		//for constant extrapolation f(i, j) = 0
		double nx = n_i(0), ny = n_i(1);
		if (dx == 0) nx = 0;
		if (dy == 0) ny = 0;
		//using the update formula for q(i, j)
		double nxny = 0;
		if ((nx!=0) || (ny!=0)){
			if((nx>0) != (ny>0)){
				nxny = nx/domain.dx - ny/domain.dy;
				dij = (nx*(dx)/domain.dx - ny*(dy)/domain.dy)/nxny;
				Pij = (nx*(Px)/domain.dx - ny*(Py)/domain.dy)/nxny;
				uij = (nx*(ux)/domain.dx - ny*(uy)/domain.dy)/nxny;
				vij = (nx*(vx)/domain.dx - ny*(vy)/domain.dy)/nxny;
			}
			else {
				nxny = nx/domain.dx + ny/domain.dy;
				dij = (nx*(dx)/domain.dx + ny*(dy)/domain.dy)/nxny;//(n_i(0)*(dx)/domain.dx + n_i(1)*(dy)/domain.dy)/nxny;
				Pij = (nx*(Px)/domain.dx + ny*(Py)/domain.dy)/nxny;
				uij = (nx*(ux)/domain.dx + ny*(uy)/domain.dy)/nxny;
				vij = (nx*(vx)/domain.dx + ny*(vy)/domain.dy)/nxny;
			}
		}

		//std::cout << dij << '\t' << i-2 << ", " << j-2 << '\t' << '\t' << dx << '\t' << dy << '\t' << '\t' << nxny << '\t' << nx << '\t' << ny << std::endl; //<< Pij << '\t' << uij << '\t' << vij << std::endl;

		primitive(i, j) = vector4(dij, uij, vij, Pij);
	};

	//Convert the conserved variables into primitive form
	vecarray primitive(domain.Nx+4, domain.Ny+4);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			primitive(i+2, j+2) = rb.fluid.state_function->primitiveVar(rb.fluid.U(i+2, j+2));
		}
	}

	boundary_conditions(primitive, domain);

	//Gauss Seidel iteration for alternate directions.
	//A total of four sweeps is performed over the entire computational domain
	//1) i = 1:I, j = 1:J
	//2) i = I:1, j = 1:J
	//3) i = I:1, J = J:1
	//4) i = 1:I, j = J:1

	// a narrowband of abs(phi) < cellwidth is imposed
	auto gauss_seidel = [heaviside, ghost_fluid, domain, ls, gr, normal](vecarray &primitive){
		//1) i = 1:I, j = 1:J
		for (int i=0; i<domain.Nx; i++){
			for (int j=0; j<domain.Ny; j++){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (gr.ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//2) i = I:1, j = 1:J
		for (int i=domain.Nx-1; i>=0; i--){
			for (int j=0; j<domain.Ny; j++){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (gr.ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//3) i = I:1, J = J:1
		for (int i=domain.Nx-1; i>=0; i--){
			for (int j=domain.Ny-1; j>=0; j--){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (gr.ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//1) i = 1:I, j = 1:J
		for (int i=0; i<domain.Nx; i++){
			for (int j=domain.Ny-1; j>=0; j--){
				if (heaviside(i+domain.buffer, j+domain.buffer)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (gr.ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}
	};

	gauss_seidel(primitive);

//////Implement narrow band here
	//reflect the ghostfluid states
	for (int i=0; i<domain.Nx; i++){
		for(int j=0; j<domain.Ny; j++){
			if (ls.phi(i+domain.buffer, j+domain.buffer) < 0){//&& ls.phi(i, j) > -(domain.dx + domain.dy)) //&& ls.phi(i+1, j+1) >= -(Test.domain.dx))
				vector2 velocity = Particle::velocity(domain.X(i, j), gr);
				//std::cout << velocity << std::endl;
				reflected_state(rb, primitive, i+2, j+2, normal(i, j), velocity);
			}
		}
	}

}

//-------------------
//Collision
//-------------------
//INCOMPLETE//
double norm(const vector2& v){
		return sqrt(pow(v(0), 2) + pow(v(1), 2));
}

double RigidBodies::compute_normal_force(const Particle& gr_i, double dist, const vector2& v_n){ //v_a is the relative velocity of the particles in contact
	//Kawamoto, Level Set DEM for 3d computations
	//normal force contribution from contact of node m_a on grain i
		//to add a damping factor to account for dissipation and increase timestep, 
		//fn =kδ+γ0vn 
	//a damping coefficient must needs be provided, forcontribution = - damping_coefficient*v_n; 
	double Fn_i = -gr_i.k_n*dist + gr_i.damping_coefficient*norm(v_n); 
	return Fn_i;
}

vector2 RigidBodies::compute_tangential_force(int j, Particle& gr_i, double Fn_i, double dist, double dt, const vector2& v_t, const vector2& normal){
//particle i --> gr_i's tangential force history will be updated
//particle j is the index of the particle undergoing collision with target i

	//if contact is active, we need to project the tangential spring onto the new tangential plane,
	//as the frame of reference will have shifted from the previous timestep
	vector2 spring_old = gr_i.springs.find(j)->second;
	vector2 spring_new = spring_old - spring_old.dot(normal)*normal;
	//We then need to enforce that |spring_new| = |spring_old|
	spring_new = spring_new/norm(spring_new)*norm(spring_old);

	//Alternatively, a better approach is to rotate the old spring into the new plane
		//...to be implemented

	//Next, we need to compute a test force to determine if there is sliding
	vector2 test_force_t = -gr_i.k_s*spring_new - gr_i.damping_coefficient*v_t;

	//By Coulomb's law, the static friction is related to normal force by
	double force_c_s = gr_i.miu*Fn_i;
	if (force_c_s < 0) {
		throw "Contact force is less than zero, check?";
	}
	double test_force = norm(test_force_t);

	//case: static friction
	//If the tangential force does not exceed the static Coulomb friction
	if (test_force <= force_c_s){
		//The tangential spring is incremented for use in the next timestep
		gr_i.springs[j] = spring_new + v_t*dt; 
	}

	//case: dynamic friction - sliding
	else {
		//If sliding friction is active, the tangential spring must be adjusted to a length 
		//that meets Coulomb's conditions
		//This ensures that in the next iteration, the test force will be equal to the dynamic Coulomb friction
		//Here, we have made the assumption that the dynamic coefficient of friction miu_d = miu_s.
		//By Coulomb's law, dynamic friction is <= static friction
		vector2 t = test_force_t/norm(test_force_t); //tangential unit vector
		gr_i.springs[j] = -(1./gr_i.k_s)*(force_c_s*t + gr_i.damping_coefficient*v_t);
	}

	//By limiting the spring length to be consistent with Coulomb's law, we effectively
	//cap the shear force to be a fraction of the normal force
	//in short, force_t = +min(f_c, |test_force|).

	return test_force_t; //Accumulated tangential force at current timestep
}
/*
vector2 RigidBodies::compute_tangential_force(Particle& gr_i, const vector2& v_t, double Fn_i, double dist, double dt){
	//particle i --> gr_i's tangential force history will be updated
	//dependant on the normal force
	//the incremental shear displacement s = vt*t
	//tanmgential velocity vt is given by vt = v - vn*normal i.e. total velocity - normal velocity
	//where vn = normal dot v //therefore vt = v_a - (va dot normal)normal
	//incremental shear displacement is given by relative tangential velocity*time

	vector2 shear_rate = (v_a - v_a.dot(normal)*normal);
	vector2 dFs = gr_i.k_s*(shear_rate*dt + gr_i.damping_coefficient*shear_rate); 
	//shear force on grain i due to node m_a is then updated using
	//gr_i.force_t = F_s; //When accumulating, must Fs be reset everytime a new contact occurs?

	double dFs_length = norm(dFs);
	vector2 F_s = dFs/dFs_length*fmin(dFs_length, gr_i.miu*norm(Fn_i));

	//since the tangential force is limited by coulomb's law, where the shear force is capped by the normal force (static friction),
	//vector2 Fs_i = gr_i.force_t/sqrt(gr_i.force_t(0)*gr_i.force_t(0) + gr_i.force_t(1)*gr_i.force_t(1))*fmin(sqrt(gr_i.force_t(0)*gr_i.force_t(0) + gr_i.force_t(1)*gr_i.force_t(1)), gr_i.miu*sqrt(Fn_i(0)*Fn_i(0) + Fn_i(1)*Fn_i(1)));
	
	return F_s;
}
*/

double RigidBodies::compute_torque(const Particle& gr_i, const vector2& node, const vector2& Fn_i, const vector2& Fs_i){
	//where node is the position vector of the seed on the level set boundary

	//moment contribution on grain i and grain j respectively
	//in 2d, this moment contribution lies in the z axis 
	double Mn_i = Particle::cross(node-gr_i.centre, Fn_i);
	//The moment contributed by the shear force of node m_a onto grain i is
	double Ms_i = Particle::cross(node-gr_i.centre, Fs_i);

	//The total torque is the sum of tangential and normal moments
	return Mn_i + Ms_i;
}

void RigidBodies::particle_collision(int j, Particle& gr_i, Particle& gr_j, const vector2& node, double dist, const vector2& normal, double dt){
	//relative velocity v_a of node m_a to grain j is
	vector2 v_a = gr_i.vc + Particle::cross(gr_i.w, node-gr_i.centre) - gr_j.vc - Particle::cross(gr_j.w, node-gr_j.centre);
	vector2 v_n = v_a.dot(normal)*normal; //normal must be a unit vector
	vector2 v_t = v_a - v_n;

	double force_n = compute_normal_force(gr_i, dist, v_n);
	vector2 Fn_i = force_n*normal;
	vector2 Fs_i = compute_tangential_force(j, gr_i, force_n, dist, dt, v_t, normal);
	gr_i.force = Fn_i + Fs_i;
	gr_i.torque = compute_torque(gr_i, node, Fn_i, Fs_i);
}

//if position.x > domain.Lx or < 0, ls = position.x - domain.lx
void RigidBodies::wall_collision(int k, Particle& gr_i, const vector2& node, double dist, const vector2& normal, double dt){
//Using the same linear elastic model as with particle collisions
//"penetration" distance and unit normal are specified
//////
	//int k is an index not shared by any other particles.
	//For this reason, negative integers are used to map k
	//The index for wall collisions are:
		//North = -1
		//East  = -2
		//South = -3
		//West  = -4

	//relative velocity v_a of node m_a to the wall is,
	//taking the wall as a stationary object:
	vector2 v_a = gr_i.vc + Particle::cross(gr_i.w, node-gr_i.centre);
	vector2 v_n = v_a.dot(normal)*normal; //normal must be a unit vector
	vector2 v_t = v_a - v_n;

	double force_n = compute_normal_force(gr_i, dist, v_n);
	vector2 Fn_i = force_n*normal;
	vector2 Fs_i = compute_tangential_force(k, gr_i, force_n, dist, dt, v_t, normal);
	gr_i.force = Fn_i + Fs_i;
	gr_i.torque = compute_torque(gr_i, node, Fn_i, Fs_i);
}

void RigidBodies::contact_detection(){
/////WORKING ON THIS//////
	//std::unordered_map<int, int>;
}

/*
void RigidBodies::contact_detection(const Domain2D& domain, std::vector<Particle>& particles, const std::vector<LevelSet>& levelset_collection, double dt){
// TO BE IMPROVED ON, made parallel
	int number_of_particles = static_cast<int>(particles.size());
	//for each particle in the system
	for (int i=0; i<number_of_particles; i++){ //master grain
	//Loop through every other particle
	//For each node of the particle
		for (int a=0; a<static_cast<int>(particles[i].nodes.size()); a++){
			Coordinates node_a(particles[i].nodes[a]);
		//remember to rotate/translate the nodes in the main solver
		//interpolate coordinates of nodes with slave levelset
		//if this value is < 0, contact is determined.
		//return a set of contact values // OR // calculate contact forces here?
			for (int j=0; j<number_of_particles; j++){ //slave grain
			//check for contact between grain i and slaves j
				if (i!=j){
					double dist = LevelSetMethods::interpolation_value(levelset_collection[j], domain, node_a);
					if (dist < 0){
						vector2 grad_phi = LevelSetMethods::interpolation_gradient(levelset_collection[j], domain, node_a);
						vector2 normal = grad_phi/(sqrt(grad_phi(0)*grad_phi(0) + grad_phi(1)*grad_phi(1)));
						particle_collision(particles[i], particles[j], particles[i].nodes[a], dist, normal, dt);
					}
				}
			}
*/
		//after calculating position future position, check that it is not at the boundaries.
		//if it crosses boundary, calvulate the extent of overlap and reflect the normal velocity?
		//check if any nodes exceed domain boundary. Assumes domain origin at 0, 0
		/*
			if ((node_a.x < 0) != (node_a.x > domain.Lx)){
				//std::cout << "contact with wall" << std::endl;
				//if node is less than neither less than 0 nor greater than Lx, it will return false
				vector2 normal(1, 0);
				//since dist must always be negative,
				double dist = (node_a.x < 0)*(node_a.x) + (node_a.x > domain.Lx)*(node_a.x - domain.Lx);
					//On the right boundary, the normal points in the wrong direction but distance is positive	
				wall_collision(particles[i], particles[i].nodes[a], dist, normal, dt);
			}
			//vertical boundaries
			if ((node_a.y < 0) != (node_a.y > domain.Ly)){
				//std::cout << "contact with wall" << std::endl;
				//if node is less than neither less than 0 nor greater than Lx, it will return false
				vector2 normal(0, 1);
				//since dist must always be negative,
				double dist = (node_a.y < 0)*(node_a.y) + (node_a.y > domain.Ly)*(node_a.y - domain.Ly);
					//On the right boundary, the normal points in the wrong direction but distance is positive	
				wall_collision(particles[i], particles[i].nodes[a], dist, normal, dt);
			}
		*/
			//std::cout << "node = " << node_a.x << '\t' << node_a.y << std::endl;
/*			if (node_a.x < 0){
				vector2 normal(1, 0);
				double dist = node_a.x;
				//std::cout << "dist = " << dist << std::endl;
				wall_collision(particles[i], particles[i].nodes[a], dist, normal, dt);
			}

			else if (node_a.x > domain.Lx){
				vector2 normal(-1, 0);
				double dist = domain.Lx - node_a.x;	
				//std::cout << "dist = " << dist << std::endl;
				wall_collision(particles[i], particles[i].nodes[a], dist, normal, dt);
			}

			if (node_a.y < 0){
				vector2 normal(0, 1);
				double dist = node_a.y;
				//std::cout << "dist = " << dist << std::endl;
				wall_collision(particles[i], particles[i].nodes[a], dist, normal, dt);
			}

			else if (node_a.y > domain.Ly){
				vector2 normal(0, -1);
				double dist = domain.Ly - node_a.y;	
				//std::cout << "dist = " << dist << std::endl;
				wall_collision(particles[i], particles[i].nodes[a], dist, normal, dt);
			}

		}
	}
	//std::cout << particles[0].centre.transpose() << '\t' << particles[1].centre.transpose() << std::endl;
	//note that the first optimisation is to prevent double checking of collisions
}*/

void RigidBodies::newton_euler(const LevelSet& ls, Particle& gr, const Domain2D& domain, const vector2& force, double torque, double dt){
//updates the velocities which are passed as arguments

//translational velocity, angular velocity, start and stop time. 
//the velocities provided need to be 'initial' values at tstart.
	//const vector2 force = LevelSetMethods::force(var, gr.ls, domain);
	//const double torque = LevelSetMethods::torque(var, gr.ls, domain, gr.centroid);

	double m = Particle::mass(gr, domain); //mass
	double J = Particle::moment_of_inertia(ls, gr, domain); //moment of inertia

	double C1 = gr.damping_coefficient*dt/2.;
	double C2 = 1/(1 + C1);

	//gr.vc contains the translational velocity of the previous timestep - n-1/2
	//double v_x_old = gr.vc(0);
	//double v_y_old = gr.vc(1);
	double v_x_new = C2*((1 - C1)*gr.vc(0) + dt/m*force);
	double v_y_new = C2*((1 - C1)*gr.vc(1) + dt/m*force);
	double w_new = C2*((1 - C1)*gr.w + dt/J*torque);

	gr.vc = vector2(v_x_new, v_y_new);
	gr.w = w_new;
}

void RigidBodies::update_displacements(Particle& gr, double dt){
//Updates the displacements and rotation for the next timestep
	gr.s += dt*gr.vc*dt; //total particle displacement
	gr.theta += dt*gr.w*dt; //total rotation
}

void RigidBodies::initial_conditions(demTests& Test){
	//Setting the initial fluid conditions
	Test.var.fluid.U.resize(Test.domain.Nx+4, Test.domain.Ny+4);
	Test.var.fluid.F.resize(Test.domain.Nx+2, Test.domain.Ny+2);
	Test.var.fluid.G.resize(Test.domain.Nx+2, Test.domain.Ny+2);


//////////////needs optimising/////////////
	Test.var.combinedls = Particle::merge(Test.var.particles, Test.domain);
	for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			for (int a=0; a<static_cast<int>(Test.var.particles.size()); a++){ 
				if (Test.var.particles[a].ls.phi(i+Test.domain.buffer, j+Test.domain.buffer) >= 0){
					if (Test.interface(i, j) == false){
						Test.var.fluid.U(i+2, j+2) = Test.var.fluid.state_function->conservedVar2Dx(Test.initialL);
					}
					else {
						Test.var.fluid.U(i+2, j+2) = Test.var.fluid.state_function->conservedVar2Dx(Test.initialR);
					}
				}

				else {
					Test.var.fluid.U(i+2, j+2).setZero();
				}
			}
			Test.var.fluid.F(i+1, j+1).setZero();
			Test.var.fluid.G(i+1, j+1).setZero();
		}
	}

	MUSCL::boundary_conditions_reflective(Test.var.fluid, Test.domain);

	//compute the normal vector using the levelset function
	std::vector<Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> > normal;
	for (int a=0; a<static_cast<int>(Test.var.particles.size()); a++){
		Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> n(Test.domain.Nx, Test.domain.Ny);
		for (int i=0; i<Test.domain.Nx; i++){
			for (int j=0; j<Test.domain.Ny; j++){
				n(i, j) = LevelSetMethods::normal(Test.var.particles[a].ls, Test.domain, i+Test.domain.buffer, j+Test.domain.buffer);
			}
		}
		normal.push_back(n);
	}

	//populate the rigid body with ghost values
	//note that these values need to be reflected before use
	//values are currently reflected in the fast sweep function
	//note: because each particle has its own unique velocity, the fast sweeping needs to be done independantly?
	for (int a=0; a<static_cast<int>(Test.var.particles.size()); a++){
		fast_sweep(Test.var.particles[a].ls, Test.var.particles[a], Test.var, Test.domain, normal[a]);
	}

	Test.var.combinedls = Particle::merge(Test.var.particles, Test.domain);
}

//can subcycling work? what of the fluid forces on particles during the contact phase?
void RigidBodies::subcycling(Moving_RB &var, const Domain2D& domain, std::vector<LevelSet>& levelset_collection, double tstop, double subdt){
	//t refers to the current timestep, tstop the timestep at the end of subcycling
	//tstop is the current timestep
	//start subcycling from the previous timestep t - dt
	double t = tstop - domain.dt;
	do{
		if (t + subdt > tstop) subdt = tstop - t;
		t += subdt;
		//update force and torque of each particle
		RigidBodies::contact_detection(domain, var.particles, levelset_collection, subdt);

		//Forces from fluids
		//For each particle:
		//forces from fluid pressure can be incorporated here. calculate before and treat as constant?
		//or calculate as particle moves along fluid, without moving the fluid?
		for (int a=0; a<static_cast<int>(var.particles.size()); a++){
			//calculate forces from fluid pressure
			//vector2 force = LevelSetMethods::force(var.fluid, levelset_collection[a], domain);
			//double torque = LevelSetMethods::torque (var.fluid, levelset_collection[a], domain, var.particles[a].centre);
			//if (count%10==0) {std::cout << force.transpose() << '\t' << torque << std::endl;}
			//calculate and update particle velocities
			newton_euler(levelset_collection[a], var.particles[a], domain, var.particles[a].force, var.particles[a].torque, domain.dt);
			update_displacements(var.particles[a], domain.dt);
			//if (count%10==0) std::cout << "particle " << a << '\t' << var.particles[a].vc.transpose() << '\t' << var.particles[a].w << std::endl;
			//if (count%10==0) std::cout << var.particles[a].force.transpose() << '\t' << var.particles[a].torque << std::endl;
//////////////ERROR///
			//move the particles
			levelset_collection[a] = var.particles[a].motion(domain);
			//update particle with new centre of gravity
			var.particles[a].centre = Particle::center_of_mass(levelset_collection[a], var.particles[a], domain);
			//extrapolate ghost values
				//compute the normal vector using the levelset function
				Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal(domain.Nx, domain.Ny);
				for (int i=0; i<domain.Nx; i++){
					for (int j=0; j<domain.Ny; j++){
						normal(i, j) = LevelSetMethods::normal(levelset_collection[a], domain, i+1, j+1);
					}
				}
			//update the ghost values in new level set position
			fast_sweep(levelset_collection[a], var.particles[a], var, domain, normal);
		}

		//std::cout << var.particles[0].force.transpose() << '\t' << var.particles[1].force.transpose() << std::endl;
	}while(t < tstop);
}

void RigidBodies::solver(Moving_RB &var, Domain2D &domain, double CFL){
	//In addition to the fluid solver of the stationary rigid body problem,
	//Forces on the rigid body are to be calculated at everytime step
	//This is an ODE problem which is solved to obtain the velocity (translational and rotational) of the rigid body
	//These velocities are then used to reimplement the level set at the new location
	//At the next time-step, boundary conditions are computed for the fluid with the new velocity and position of the rigid body
	
	//Storage for (piecewise linear) interpolated cell values 
	matrix ULix(domain.Nx+4, 4);
	matrix URix(domain.Nx+4, 4);

	matrix ULiy(domain.Ny+4, 4);
	matrix URiy(domain.Ny+4, 4);

	//Temporary storage for x and y conserved variables 
	matrix Uxn(domain.Nx+4, 4);
	matrix Uyn(domain.Ny+4, 4);


	slopeLimiter a = VanLeer;//getLimiter();
	/////-----
	//LevelSets//
	//---------------------
	//When solving for the fluid values, we can consider all the levelsets at once by taking their union
	var.combinedls = Particle::merge(var.particles, domain);
	//storing the initial collection of levelsets
	std::vector<LevelSet> levelset_collection;
	for (int a=0; a<static_cast<int>(var.particles.size()); a++){
		levelset_collection.push_back(var.particles[a].ls);
	}

	double faket = 0.1;
	double t = 0.0;
	double plot_time = 0.1; //0.05; //0.1;
	int count = 0;
	do{
		//set timestep, taking maximum wavespeed in x or y directions
		//looping through the entire domain, find the maximum possible wavespeed in both directions
		double Smax=0;
		//double Smax_y=0;
		for (int i=2; i<domain.Nx+2; i++){
			for (int j=2; j<domain.Ny+2; j++){
				//checks if density is broken
				if (var.fluid.U(i, j)(0) != var.fluid.U(i, j)(0) || var.fluid.U(i, j)(0) < 0){
					std::cout << "Density value is not physical - NaN, terminating" << std::endl;
					throw "Density value is not physical";
				}
				vector4 Utmp = var.fluid.U(i, j);
				double a_ij = var.fluid.state_function->soundspeed(Utmp);
				double d_ij_x = Utmp(1)/Utmp(0); //velocity in x direction
				double d_ij_y = Utmp(3)/Utmp(0); //velocity in y direction
				double vn = sqrt(pow(d_ij_x, 2) + pow(d_ij_y, 2));
				double lambda = fmax(fmax(abs(d_ij_x) + a_ij, abs(d_ij_y) + a_ij), vn + a_ij);
				if (lambda > Smax) Smax = lambda;
			}
		}

		if (count <= 5){
				//using a CFL number of 0,.2 for the first 5 timesteps
				domain.dt = 0.2*fmin((domain.dx/Smax), (domain.dy/Smax));
			}
		else {
				domain.dt = CFL*fmin((domain.dx/Smax), (domain.dy/Smax));
		}
		if (t + domain.dt > plot_time) {
			domain.dt = plot_time - t;
			t = plot_time;
			//std::cout << t << "s plotting" << std::endl;
		}
		//if (t + domain.dt > domain.tstop) domain.dt = domain.tstop - t;
		else{
			t += domain.dt;
		}
		count += 1;

		//sweeping in the x-direction for each y row
		for (int j=0; j<domain.Ny; j++){
			//Within each row, store the values for each grid point in x with a single row list
			for(int i=0; i<domain.Nx+4; i++){
				Uxn.row(i) = var.fluid.U(i, j+2);
				//the sweep in x is only performed within the real y grid
			} 
			//calculate and update to new interpolated values U_i
			MUSCL::data_reconstruction(Uxn, a, ULix, URix, domain.Nx);

			//sweeping in the x-direction for each y row
			for (int i=1; i<domain.Nx+2; i++){
				if (var.combinedls.phi(i+domain.buffer-1, j+domain.buffer) >= 0 || var.combinedls.phi(i+domain.buffer, j+domain.buffer) >= 0 || var.combinedls.phi(i+domain.buffer-2, j+domain.buffer) >= 0){
					vector4 Fx(0, 0, 0, 0);
					//edge case where there is only one ghost cell
					if (var.combinedls.phi(i+domain.buffer, j+domain.buffer) < 0 && var.combinedls.phi(i+domain.buffer+1, j+domain.buffer) >= 0){
						MUSCL::compute_fluxes(var.fluid, Uxn, Fx, i);
						//std::cout << i-1 << '\t' << j << '\t' << var.levelsets[0].phi(i+1, j+1) << std::endl;
					}
					//normal case
					else {
						MUSCL::compute_fluxes(var.fluid, domain, Uxn, Fx, i, ULix, URix, domain.dx);
						//std::cout << Fx.transpose() << '\t' << '(' << i << ", " << j+1 << ')' << std::endl;
						//if (i==13) std::cout << Uxn.row(i) << '\t' << '\t' << ULix.row(i) << '\t' << '\t'  << URix.row(i) << std::endl;
					}
					var.fluid.F(i, j+1) = Fx; //storing the computed flux
				}
			}
		}

	//for (int i=0; i<domain.Nx; i++){
	//	for (int j=0; j<domain.Ny; j++){
			//std::cout << var.fluid.F(i+1, j+1).transpose() << std::endl;
	//	}
	//}

		//sweeping in the y-direction for each x row
		for (int i=0; i<domain.Nx; i++){
			for (int j=0; j<domain.Ny+4; j++){
				//Since the normal velocity is now v, its position is swapped with u
				Uyn.row(j) = var.fluid.swap_xy(var.fluid.U(i+2, j));
			}
			MUSCL::data_reconstruction(Uyn, a, ULiy, URiy, domain.Ny);

			for (int j=1; j<domain.Ny+2; j++){
				if (var.combinedls.phi(i+domain.buffer, j+domain.buffer-1) >= 0 || var.combinedls.phi(i+domain.buffer, j+domain.buffer+1) >= 0 || var.combinedls.phi(i+domain.buffer, j+domain.buffer-2) >=0){
					vector4 Fy(0, 0, 0, 0);
					if (var.combinedls.phi(i+domain.buffer, j+domain.buffer) < 0 && var.combinedls.phi(i+domain.buffer, j+domain.buffer+1) >= 0){
						MUSCL::compute_fluxes(var.fluid, Uyn, Fy, j);
						//std::cout << i << '\t' << j-1 << '\t' << var.levelsets[0].phi(i+1, j+1) << std::endl;
					}	
					else {
						MUSCL::compute_fluxes(var.fluid, domain, Uyn, Fy, j, ULiy, URiy, domain.dy);
					}
					var.fluid.G(i+1, j) = Fy; //storing the computed flux.
				}		
			}
		}
/*
		if (count >= 100){
			for (int i=0; i<domain.Nx; i++){
				for (int j=0; j<domain.Nx; j++){
					//std::cout << var.fluid.F(i+1, j+1).transpose() << '\t' << '(' << i << ", " << j << ')' << std::endl;
					if (i==35 && j==55) {
						std::cout << "count = " << count << std::endl;
						std::cout << var.fluid.U(i+2, j+2).transpose() << " ( " << i << ", " << j << " ) " << var.combinedls.phi(i+1, j+1) << std::endl;
						std::cout << var.fluid.U(i+1, j+2).transpose() << " ( " << i-1 << ", " << j << " ) " << var.combinedls.phi(i, j+1) << std::endl;
						std::cout << var.fluid.U(i+3, j+2).transpose() << " ( " << i+1 << ", " << j << " ) " << var.combinedls.phi(i+2, j+1) << std::endl;
						std::cout << var.fluid.U(i+2, j+1).transpose() << " ( " << i << ", " << j-1 << " ) " << var.combinedls.phi(i+1, j) << std::endl;
						std::cout << var.fluid.U(i+2, j+3).transpose() << " ( " << i << ", " << j+1 << " ) " << var.combinedls.phi(i+1, j+2) << std::endl;
						std::cout << var.fluid.F(i+1, j+1).transpose() << '\t' << '(' << i << ", " << j << ')' << std::endl;
					}
					//std::cout << var.fluid.state_function->primitiveVar(var.fluid.U(i+2, j+2)).transpose() << " ( " << i << ", " << j << " ) " << std::endl;
					//if (var.fluid.state_function->primitiveVar(var.fluid.U(i+2, j+2))(3) < 0 || var.fluid.state_function->primitiveVar(var.fluid.U(i+2, j+2))(3) > 5) {
					//	std::cout << var.fluid.state_function->primitiveVar(var.fluid.U(i+2, j+2)).transpose() << " ( " << i << ", " << j << " ) " << std::endl;
					//}
				}
			}
		}
*/

		//updating U with Strand splitting -- X(0.5dt)Y(dt)X(0.5dt)
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				if (var.combinedls.phi(i+domain.buffer-2, j+domain.buffer) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i, j+2), var.fluid.F(i, j+1), var.fluid.F(i-1, j+1), domain.dt/2, domain.dx);
					//std::cout << var.fluid.F(i, j+1).transpose()  << '\t'  << '\t' << var.fluid.F(i-1, j+1).transpose() << std::endl;
				}
			}
		}
		for (int i=0; i<domain.Nx; i++){
			for (int j=2; j<domain.Ny+2; j++){
				if (var.combinedls.phi(i+domain.buffer, j+domain.buffer-2) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i+2, j), var.fluid.swap_xy(var.fluid.G(i+1, j)), var.fluid.swap_xy(var.fluid.G(i+1, j-1)), domain.dt, domain.dy);
				}
			}
		}		
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				if (var.combinedls.phi(i+domain.buffer-2, j+domain.buffer) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i, j+2), var.fluid.F(i, j+1), var.fluid.F(i-1, j+1), domain.dt/2, domain.dx);
				}
			}
		}
		MUSCL::boundary_conditions_reflective(var.fluid, domain);

		//---------------------------------------------
		//	Forces and Torque on the Particle, DEM
		//---------------------------------------------
		RigidBodies::subcycling(var, domain, levelset_collection, t, 0.01);
/*
		//Perform contact detection, accumulating contact forces on the particles
		//RigidBodies::contact_detection(domain, var.particles, levelset_collection, domain.dt);
		//reset the levelset collection 
		
		//Forces from fluids
		for (int a=0; a<static_cast<int>(var.particles.size()); a++){
			//calculate forces from fluid pressure
			//vector2 force = LevelSetMethods::force(var.fluid, levelset_collection[a], domain);
			//double torque = LevelSetMethods::torque (var.fluid, levelset_collection[a], domain, var.particles[a].centre);
			//if (count%10==0) {std::cout << force.transpose() << '\t' << torque << std::endl;}
			//calculate and update particle velocities
//////////////ERROR///
//unstable unless low cfl is used.. tiomestep issues
			newton_euler(levelset_collection[a], var.particles[a], domain, var.particles[a].force, var.particles[a].torque, domain.dt);
			//if (count%10==0) std::cout << "particle " << a << '\t' << var.particles[a].vc.transpose() << '\t' << var.particles[a].w << std::endl;
			//if (count%10==0) std::cout << var.particles[a].force.transpose() << '\t' << var.particles[a].torque << std::endl;
//////////////ERROR///
			//move the particles
			levelset_collection[a] = var.particles[a].motion(domain, t);
			//update particle with new centre of gravity
			var.particles[a].centre = Particle::center_of_mass(levelset_collection[a], var.particles[a], domain);
			//extrapolate ghost values
				//compute the normal vector using the levelset function
				Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal(domain.Nx, domain.Ny);
				for (int i=0; i<domain.Nx; i++){
					for (int j=0; j<domain.Ny; j++){
						normal(i, j) = LevelSetMethods::normal(levelset_collection[a], domain, i+1, j+1);
					}
				}
			fast_sweep(levelset_collection[a], var.particles[a], var, domain, normal);
		}
*/

		//Merge the new level sets
		var.combinedls = LevelSetMethods::merge(levelset_collection, domain);
		
		if (count%50==0) std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << domain.dt << std::endl;
		//if (count == 60) t = domain.tstop;
		//if (count == 10) std::cout << t << std::endl;
		if (t == plot_time) {
			std::cout << "plotting " << t << std::endl;
			std::string file = "Data/dataeuler_" + std::to_string(faket) + ".txt";
			std::string file2 = "Data/datapoints_" + std::to_string(faket) + ".txt";
			output(var, domain, file, file2);
			std::cout << "particle " << 0 << '\t' << var.particles[0].vc.transpose() << '\t' << var.particles[0].w << std::endl;
			std::cout << "particle " << 1 << '\t' << var.particles[1].vc.transpose() << '\t' << var.particles[1].w << std::endl;
			plot_time += 0.1;//0.1;
			faket+=0.1;
		}

	}while (t < domain.tstop);
	std::cout << "count = "  << count << std::endl;
}

void RigidBodies::output(const Moving_RB &var, const Domain2D &domain, std::string filename, std::string filename2){

	std::ofstream outfile;
	std::ofstream outfile2;

	outfile.open(filename);
	outfile2.open(filename2);

	for (int z=0; z<static_cast<int>(var.particles[0].nodes.size()); z++){
		//std::cout << var.particles[0].nodes[z].transpose() << std::endl;
		outfile2 << var.particles[0].nodes[z].transpose() << std::endl<< std::endl;
	}

	double u=0;
	double P=0;
	double e=0;
	double schlieren=0;

	for (int i=2; i<domain.Nx+2; i++){
		for (int j=2; j<domain.Ny+2; j++){
			if (var.combinedls.phi(i+domain.buffer-2, j+domain.buffer-2) >= 0){
				vector4 Ux = var.fluid.U(i, j);
				if (Ux(0)!=0){
					u = sqrt(pow(Ux(1)/Ux(0), 2) + pow(Ux(3)/Ux(0), 2));
					P = var.fluid.state_function->Pressure(Ux);
					e = var.fluid.state_function->internalE(Ux);

					//central difference to calculate partial derivatives in x, y for density
					double grad_density_x = (var.fluid.U(i+1, j)(0) - var.fluid.U(i-1, j)(0))/(2*domain.dx);
					double grad_density_y = (var.fluid.U(i, j+1)(0) - var.fluid.U(i, j-1)(0))/(2*domain.dy);
					//calculating the numerical schlieren
					schlieren = exp((-20*sqrt(pow(grad_density_x, 2) + pow(grad_density_y, 2)))/(1000*Ux(0)));
				}

				outfile << domain.dx*(i-2) << '\t' << domain.dy*(j-2) << '\t' << Ux(0) << '\t' << u
				<< '\t' << P << '\t' << e << '\t' << schlieren << '\t' << '\t' << '\t' << var.combinedls.phi(i+domain.buffer-2, j+domain.buffer-2) << '\t' << std::endl;
			}
			/*else {
				outfile << domain.dx*(i-2) << '\t' << domain.dy*(j-2) << '\t' << 0 << '\t' << 0
				<< '\t' << 0 << '\t' << 0 << '\t' << 0 << std::endl;
			}*/
		}
		outfile << std::endl;
	}

	outfile.close();
	outfile2.close();
}

void RigidBodies::output_levelset(const Moving_RB &var, const Domain2D &domain, std::string filename, std::string filename2){

	std::ofstream outfile;
	std::ofstream outfile2;

	outfile.open(filename);
	outfile2.open(filename2);

	for (int z=0; z<static_cast<int>(var.particles[0].nodes.size()); z++){
		//std::cout << var.particles[0].nodes[z].transpose() << std::endl;
		outfile2 << var.particles[0].nodes[z].transpose() << std::endl<< std::endl;
	}


	for (int i=2; i<domain.Nx+2; i++){
		for (int j=2; j<domain.Ny+2; j++){
			outfile << domain.dx*(i-2) << '\t' << domain.dy*(j-2) << '\t' << var.combinedls.phi(i-1, j-1) << '\t' << std::endl;
		}
		outfile << std::endl;
	}

	outfile.close();
	outfile2.close();
}

void RigidBodies::rigid_body_solver(demTests &Test, double CFL){
	initial_conditions(Test);
	//for (int b=0; b<static_cast<int>(var.particles[a].nodes.size()); b++){
	//	std::cout << var.particles[a].nodes[b].transpose() << std::endl;
	//}
	//std::cout << Test.var.particles[0].nodes.size() << std::endl;
	solver(Test.var, Test.domain, CFL);
	output(Test.var, Test.domain, "dataeuler.txt", "datapoints.txt");
	std::cout << "done: Rigid Body" << std::endl;
}





