#include "../headerFiles/RigidBodies.h"

void RigidBodies::fast_sweep(const LevelSet &ls, RB_2D &rb, const Domain2D &domain, const Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> &normal){

	//Multi dimensional extrapolation using fast sweeping
	auto heaviside = [ls](int i, int j){
		if (ls.phi(i, j) >= 0) return 0;
		else return 1;
	};

	auto ghost_fluid = [rb, domain](int i, int j, vector2 n_i){
		//info is swept from outside into the rigid body
		//constant extrapolation of density (d), pressure (P) and velocity (v_i)
		//using the first order upwind finite difference approximation
		double dij=0, Pij=0, uij=0, vij=0;
		double dx, Px=0, ux=0, vx=0;
		double dy, Py=0, uy=0, vy=0;
		//consider the x derivatives
		if (n_i(0) > 0){
			//information ravels from positive to negative phi.
			//if normal > 0, 
			//information is travelling to the left
			//right side is upwind
			dx = rb.fluid.U(i+1,j)(0);
			if(dx>0){
				Px = rb.fluid.state_function->Pressure(rb.fluid.U(i+1,j));
				ux = rb.fluid.U(i+1,j)(1)/rb.fluid.U(i+1,j)(0); //u
				vx = rb.fluid.U(i+1,j)(3)/rb.fluid.U(i+1,j)(0); //v
			}
		}
		else {
			//information is travelling to the right
			//left side is upwind
			dx = rb.fluid.U(i-1,j)(0);
			if(dx>0){
				Px = rb.fluid.state_function->Pressure(rb.fluid.U(i-1,j));
				ux = rb.fluid.U(i-1,j)(1)/rb.fluid.U(i-1,j)(0);	
				vx = rb.fluid.U(i-1,j)(3)/rb.fluid.U(i-1,j)(0); 
			}
		}

		if (n_i(1) > 0){
			dy = rb.fluid.U(i,j+1)(0);
			if(dy>0){
				Py = rb.fluid.state_function->Pressure(rb.fluid.U(i,j+1));
				uy = rb.fluid.U(i,j+1)(1)/rb.fluid.U(i,j+1)(0); //u
				vy = rb.fluid.U(i,j+1)(3)/rb.fluid.U(i,j+1)(0); //v
			}
		}

		else {
			dy = rb.fluid.U(i,j-1)(0);
			if(dy>0){
				Py = rb.fluid.state_function->Pressure(rb.fluid.U(i,j-1));
				uy = rb.fluid.U(i,j-1)(1)/rb.fluid.U(i,j-1)(0); //u
				vy = rb.fluid.U(i,j-1)(3)/rb.fluid.U(i,j-1)(0); //v
			}
		}
		//for constant extrapolation f(i, j) = 0
		//using the update formula for q(i, j)
		double nxny = n_i(0)/domain.dx + n_i(1)/domain.dy;
		if (nxny !=0){
			dij = (n_i(0)*(dx)/domain.dx + n_i(1)*(dy)/domain.dy)/nxny;
			Pij = (n_i(0)*(Px)/domain.dx + n_i(1)*(Py)/domain.dy)/nxny;
			uij = (n_i(0)*(ux)/domain.dx + n_i(1)*(uy)/domain.dy)/nxny;
			vij = (n_i(0)*(vx)/domain.dx + n_i(1)*(vy)/domain.dy)/nxny;
		}

		//std::cout << nxny << '\t' << n_i(0) << '\t' << n_i(1) << std::endl;
		//std::cout << dx << '\t' << Px << '\t' << ux << '\t' << vx << std::endl;
		//std::cout << dij << '\t' << Pij << '\t' << uij << '\t' << vij << std::endl;

		vector4 extrapolatedW(dij, uij, vij, Pij);
		vector4 extrapolatedU = rb.fluid.state_function->conservedVar2Dx(extrapolatedW);
		return extrapolatedU;
	};

	//Gauss Seidel iteration for alternate directions.
	//A total of four sweeps is performed over the entire computational domain
	//1) i = 1:I, j = 1:J
	//2) i = I:1, j = 1:J
	//3) i = I:1, J = J:1
	//4) i = 1:I, j = J:1

	// a narrowband of abs(phi) < cellwidth is imposed

	//1) i = 1:I, j = 1:J
	for (int i=0; i<domain.Nx; i++){
		for (int j=0; j<domain.Ny; j++){
			if (heaviside(i+1, j+1)){
				rb.fluid.U(i+2, j+2) = ghost_fluid(i+2, j+2, normal(i, j));
				if (ls.phi(i, j) > -(domain.dx)) break;
			}
		}
	}

	//2) i = I:1, j = 1:J
	for (int i=domain.Nx-1; i>=0; i--){
		for (int j=0; j<domain.Ny; j++){
			if (heaviside(i+1, j+1)){
				rb.fluid.U(i+2, j+2) = ghost_fluid(i+2, j+2, normal(i, j));
				if (ls.phi(i, j) > -(domain.dx)) break;
			}
		}
	}

	//3) i = I:1, J = J:1
	for (int i=domain.Nx-1; i>=0; i--){
		for (int j=domain.Ny-1; j>=0; j--){
			if (heaviside(i+1, j+1)){
				rb.fluid.U(i+2, j+2) = ghost_fluid(i+2, j+2, normal(i, j));
				if (ls.phi(i, j) > -(domain.dx)) break;
			}
		}
	}

	//1) i = 1:I, j = 1:J
	for (int i=0; i<domain.Nx; i++){
		for (int j=domain.Ny-1; j>=0; j--){
			if (heaviside(i+1, j+1)){
				rb.fluid.U(i+2, j+2) = ghost_fluid(i+2, j+2, normal(i, j));
				if (ls.phi(i, j) > -(domain.dx)) break;
			}
		}
	}
}

vector4 RigidBodies::reflected_state(const RB_2D &var, int i, int j, const vector2 &n_i){
	//to be used when (phi<0) --> ghost fluid state within rigid body

	//extrapolated values
	double density = var.fluid.U(i, j)(0);
	double pressure = var.fluid.state_function->Pressure(var.fluid.U(i, j));
	vector2 velocity(var.fluid.U(i, j)(1)/var.fluid.U(i, j)(0), var.fluid.U(i, j)(3)/var.fluid.U(i, j)(0)); //(u, v)

	vector2 v_ghost = velocity - 2*(n_i.dot(velocity))*(n_i);
	//reflected state, conserved variables
	vector4 reflected(density, v_ghost(0), v_ghost(1), pressure);
	vector4 reflected_cons = var.fluid.state_function->conservedVar2Dx(reflected);

	return reflected_cons;
}

void RigidBodies::boundary_conditions(RB_2D &var, const Domain2D &domain){
	//assigning ghost values in the x-direction 
	for (int j=0; j<domain.Ny; j++){
		var.fluid.U(1, j+2) = var.fluid.U(2, j+2);
		var.fluid.U(0, j+2) = var.fluid.U(1, j+2);
		var.fluid.U(domain.Nx+2, j+2) = var.fluid.U(domain.Nx+1, j+2);
		var.fluid.U(domain.Nx+3, j+2) = var.fluid.U(domain.Nx+2, j+2);
	} 
	//assigning ghost values in the y-direction
	for (int i=0; i<domain.Nx; i++){
		var.fluid.U(i+2, 1) = var.fluid.U(i+2, 2);
		var.fluid.U(i+2, 0) = var.fluid.U(i+2, 1);
		var.fluid.U(i+2, domain.Ny+2) = var.fluid.U(i+2, domain.Ny+1);
		var.fluid.U(i+2, domain.Ny+3) = var.fluid.U(i+2, domain.Nx+2);
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
		}
	}

	boundary_conditions(Test.var, Test.domain);

	//compute the normal vector using the levelset function
	Eigen::Array<vector2, Eigen::Dynamic, Eigen::Dynamic> normal(Test.domain.Nx, Test.domain.Ny);
	for (int i=0; i<Test.domain.Nx; i++){
		for (int j=0; j<Test.domain.Ny; j++){
			normal(i, j) = LevelSetMethods::normal(Test.var.levelsets[0], Test.domain, i+1, j+1);
		}
	}
	//populate the rigid body with ghost values
	//note that these values need to be reflected before use
	for (int a=0; a<static_cast<int>(Test.var.levelsets.size()); a++){
		int count = 0;
		do{
			fast_sweep(Test.var.levelsets[a], Test.var, Test.domain, normal);
			count++;
		}while(count < 4);
	}

	std::ofstream outfile;
	outfile.open("extrapolation.txt");
	for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			outfile << i*Test.domain.dx << '\t' << j*Test.domain.dy << '\t' << Test.var.fluid.U(i+2, j+2)(0) << '\t' << Test.var.fluid.U(i+2, j+2)(1) << '\t' << Test.var.fluid.U(i+2, j+2)(2) << '\t' << Test.var.fluid.U(i+2, j+2)(3)<< std::endl;
		}
		outfile << std::endl;
	}
	outfile.close();
	
	for(int i=0; i<Test.domain.Nx; i++){
		for(int j=0; j<Test.domain.Ny; j++){
			std::cout << Test.var.fluid.U(i+2, j+2)(0) << '\t';
			//std::cout << Test.var.levelsets[0].phi(i+1, j+1) << '\t';
		}
		std::cout << std::endl;
	}
}
/*
void RigidBodies::solver(){

}

void RigidBodies::output(){

}

void RigidBodies::rigid_body_solver(){

}
*/