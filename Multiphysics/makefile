CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -g

RPsolvers/Solvers.o : RPsolvers/Solvers.C RPsolvers/Solvers.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

GhostFluid.o : GhostFluid.C GhostFluid.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

RPsolvers/LevelSet.o : RPsolvers/LevelSet.C RPsolvers/LevelSet.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/EOS.o : EOS/EOS.C EOS/EOS.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/JWL.o : EOS/JWL.C EOS/JWL.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

Tests/eulerTests.o : Tests/eulerTests.C Tests/eulerTests.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

practical2-ghostfluid.o : practical2-ghostfluid.C
	${CXX} -c $< -o $@ ${CXXFLAGS}

prac2gfm : practical2-ghostfluid.o
	g++ $^ -lm -o $@

main.o : main.C
	${CXX} -c $< -o $@ ${CXXFLAGS}

main : main.o  RPsolvers/Solvers.o GhostFluid.o RPsolvers/LevelSet.o EOS/EOS.o EOS/JWL.o Tests/eulerTests.o
	g++ $^ -lm -o $@

