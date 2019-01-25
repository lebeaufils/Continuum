CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -g

RPsolvers/Solvers.o : RPsolvers/Solvers.C RPsolvers/Solvers.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/IG.o : EOS/IG.C EOS/IG.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/JWL.o : EOS/JWL.C EOS/JWL.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

Tests/eulerTests.o : Tests/eulerTests.C Tests/eulerTests.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

main.o : main.C
	${CXX} -c $< -o $@ ${CXXFLAGS}

main : main.o  RPsolvers/Solvers.o EOS/IG.o EOS/JWL.o Tests/eulerTests.o
	g++ $^ -lm -o $@

