# Rigid-Fluid simulations

This is a numerical model that simulates the interactions between fluid and rigid bodies (particulate systems).

## Getting Started

### Prerequisites

Boost, Eigen, C++11

### Running simulations

Sadly, the code needs to be recompiled every time a parameter/test is changed.
There are 6 existing tests that can be run by changing the test case in main.C
Resolution of the domain is determined when constructing the test class
```
Test demTests(N) or Test demTests(Nx, Ny),  where N is the number of grid cells
```
Then, to run the tests, run the respective test as a function
```
Test.test1();
```
Finally, to call the solver
```
RigidBodies::rigid_body_solver(Test, CFL, plot_time, dtmax);
```
Where Test is the initialised test case and CFL is the Courant number for the fluid solver (typically 0.5).
plot_time is a double that determines the time interval data will be printed for intermediate states.
Dtmax represents the maximum time step allowed in the system, for both fluid and particles.
For the dem models to be stable, this is usually 0.002/0.001.

For example, to run the configured test1, in main.C,
```
	demTests Tests(401); 
	Tests.test1();
	RigidBodies::rigid_body_solver(Tests, 0.5, 0.1, 0.002);
```

## Examples
