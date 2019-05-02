#include "../headerFiles/RigidBodies.h"

void RigidBodies::reflected_state(RB_2D &var, const vecarray &primV, int i, int j, const vector2 &n_i){
	vector2 velocity(primV(i,j)(1), primV(i,j)(2));
	vector2 v_ghost = velocity - 2*(n_i.dot(velocity))*(n_i);

	vector4 reflected(primV(i,j)(0), v_ghost(0), v_ghost(1), primV(i,j)(3));
	vector4 reflected_cons = var.fluid.state_function->conservedVar2Dx(reflected);

	var.fluid.U(i, j) = reflected_cons;
}

void RigidBodies::fast_sweep(const LevelSet &ls, RB_2D &rb, const Domain2D &domain, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> &normal){

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
	vecarray primitive(domain.Nx+4, domain.Ny+4);
	for (int i=0; i<domain.Nx+4; i++){
		for (int j=0; j<domain.Ny+4; j++){
			primitive(i, j) = rb.fluid.state_function->primitiveVar(rb.fluid.U(i, j));
		}
	}


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
				if (heaviside(i+1, j+1)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//2) i = I:1, j = 1:J
		for (int i=domain.Nx-1; i>=0; i--){
			for (int j=0; j<domain.Ny; j++){
				if (heaviside(i+1, j+1)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//3) i = I:1, J = J:1
		for (int i=domain.Nx-1; i>=0; i--){
			for (int j=domain.Ny-1; j>=0; j--){
				if (heaviside(i+1, j+1)){
					ghost_fluid(primitive, i+2, j+2, normal(i, j));
					//if (ls.phi(i+1, j+1) < -(domain.dx + domain.dy)) break;
				}
			}
		}

		//1) i = 1:I, j = 1:J
		for (int i=0; i<domain.Nx; i++){
			for (int j=domain.Ny-1; j>=0; j--){
				if (heaviside(i+1, j+1)){
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
//////////////needs optimising/////////////
	for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			for (int a=0; a<static_cast<int>(Test.var.levelsets.size()); a++){ 
				if (Test.var.levelsets[a].phi(i+1, j+1) >= 0){
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

	MUSCL::boundary_conditions(Test.var.fluid, Test.domain);

	//compute the normal vector using the levelset function
	Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal(Test.domain.Nx, Test.domain.Ny);
	for (int i=0; i<Test.domain.Nx; i++){
		for (int j=0; j<Test.domain.Ny; j++){
			normal(i, j) = LevelSetMethods::normal(Test.var.levelsets[0], Test.domain, i+1, j+1);
		}
	}

	//populate the rigid body with ghost values
	//note that these values need to be reflected before use
	//values are currently reflected in the fast sweep function
	for (int a=0; a<static_cast<int>(Test.var.levelsets.size()); a++){
		fast_sweep(Test.var.levelsets[a], Test.var, Test.domain, normal);
	}

/*	std::ofstream outfile;
	outfile.open("extrapolation.txt");
	for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			outfile << i*Test.domain.dx << '\t' << j*Test.domain.dy << '\t' << Test.var.fluid.U(i+2, j+2)(0) << '\t' << Test.var.fluid.U(i+2, j+2)(1) << '\t' << Test.var.fluid.U(i+2, j+2)(2) << '\t' << Test.var.fluid.U(i+2, j+2)(3)<< std::endl;
		}
		outfile << std::endl;
	}
	outfile.close();
*/	
	/*for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			std::cout << Test.var.fluid.U(i+2, j+2)(2) << '\t';
			//std::cout << Test.var.levelsets[0].phi(i+1, j+1) << '\t';
		}
		std::cout << std::endl;
	}*/
}

void RigidBodies::solver(RB_2D &var, Domain2D &domain, double CFL){
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

	//compute the normal vector using the levelset function. only valid for non moving levelsets
	Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal(domain.Nx, domain.Ny);
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			normal(i, j) = LevelSetMethods::normal(var.levelsets[0], domain, i+1, j+1);
		}
	}

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
				if (var.levelsets[0].phi(i, j+1) >= 0 || var.levelsets[0].phi(i+1, j+1) >= 0 || var.levelsets[0].phi(i-1, j+1) >= 0){
					vector4 Fx(0, 0, 0, 0);
					//edge case where there is only one ghost cell
					if (var.levelsets[0].phi(i+1, j+1) < 0 && var.levelsets[0].phi(i+2, j+1) >= 0){
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


		//sweeping in the y-direction for each x row
		for (int i=0; i<domain.Nx; i++){
			for (int j=0; j<domain.Ny+4; j++){
				//Since the normal velocity is now v, its position is swapped with u
				Uyn.row(j) = var.fluid.swap_xy(var.fluid.U(i+2, j));
			}
			MUSCL::data_reconstruction(Uyn, a, ULiy, URiy, domain.Ny);

			for (int j=1; j<domain.Ny+2; j++){
				if (var.levelsets[0].phi(i+1, j) >= 0 || var.levelsets[0].phi(i+1, j+1) >= 0 || var.levelsets[0].phi(i+1, j-1) >=0){
					vector4 Fy(0, 0, 0, 0);
					if (var.levelsets[0].phi(i+1, j+1) < 0 && var.levelsets[0].phi(i+1, j+2) >= 0){
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

		//updating U with Strand splitting -- X(0.5dt)Y(dt)X(0.5dt)
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				if (var.levelsets[0].phi(i-1, j+1) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i, j+2), var.fluid.F(i, j+1), var.fluid.F(i-1, j+1), domain.dt/2, domain.dx);
					//std::cout << var.fluid.F(i, j+1).transpose()  << '\t'  << '\t' << var.fluid.F(i-1, j+1).transpose() << std::endl;
				}
			}
		}
		for (int i=0; i<domain.Nx; i++){
			for (int j=2; j<domain.Ny+2; j++){
				if (var.levelsets[0].phi(i+1, j-1) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i+2, j), var.fluid.swap_xy(var.fluid.G(i+1, j)), var.fluid.swap_xy(var.fluid.G(i+1, j-1)), domain.dt, domain.dy);
				}
			}
		}		
		for (int j=0; j<domain.Ny; j++){
			for (int i=2; i<domain.Nx+2; i++){
				if (var.levelsets[0].phi(i-1, j+1) >= 0){
					MUSCL::conservative_update_formula_2D(var.fluid.U(i, j+2), var.fluid.F(i, j+1), var.fluid.F(i-1, j+1), domain.dt/2, domain.dx);
				}
			}
		}
		MUSCL::boundary_conditions(var.fluid, domain);

		//Extrapolate new ghost values
		fast_sweep(var.levelsets[0], var, domain, normal);

		if (count%50==0) std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << "count = " << count << '\t' << t << "s" << '\t' << domain.dt << std::endl;
		//std::cout << Smax_x << '\t' << Smax_y << std::endl;
		//if (count == 1) t = domain.tstop;
	}while (t < domain.tstop);
	std::cout << "count = "  << count << std::endl;
}

void RigidBodies::output(const RB_2D &var, const Domain2D &domain){

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



