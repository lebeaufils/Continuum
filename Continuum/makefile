CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -g

RPsolvers/Solvers.o : RPsolvers/Solvers.C RPsolvers/Solvers.h 
	${CXX} -c $< -o $@ ${CXXFLAGS}

#GhostFluid.o : GhostFluid.C GhostFluid.h
#	${CXX} -c $< -o $@ ${CXXFLAGS}

LevelSet/LevelSet.o : LevelSet/LevelSet.C LevelSet/LevelSet.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/EOS.o : EOS/EOS.C EOS/EOS.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

Tests/eulerTests.o : Tests/eulerTests.C Tests/eulerTests.h 
	${CXX} -c $< -o $@ ${CXXFLAGS}

main.o : main.C
	${CXX} -c $< -o $@ ${CXXFLAGS}

#main : main.o   RPsolvers/Solvers.o GhostFluid.o LevelSet/LevelSet.o EOS/EOS.o Tests/eulerTests.o
#	g++  $^ -lm -o $@

main : main.o   RPsolvers/Solvers.o LevelSet/LevelSet.o EOS/EOS.o Tests/eulerTests.o
	g++  $^ -lm -o $@

memcheck : main.o   RPsolvers/Solvers.o LevelSet/LevelSet.o EOS/EOS.o Tests/eulerTests.o
	g++ -g  $^ -o $@



