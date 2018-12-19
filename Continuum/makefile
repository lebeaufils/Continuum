CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -g

RPsolvers/HLLC.o : RPsolvers/HLLC.C RPsolvers/HLLC.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

RPsolvers/MUSCL.o : RPsolvers/MUSCL.C RPsolvers/MUSCL.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/IG.o : EOS/IG.C EOS/IG.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

EOS/JWL.o : EOS/JWL.C EOS/JWL.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

Tests/eulerTests.o : Tests/eulerTests.C Tests/eulerTests.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

main.o : main.C
	${CXX} -c $< -o $@ ${CXXFLAGS}

main : main.o RPsolvers/HLLC.o RPsolvers/MUSCL.o EOS/IG.o EOS/JWL.o Tests/eulerTests.o
	g++ $^ -lm -o $@

