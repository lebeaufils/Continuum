CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -g

objectFiles/Solvers.o : cppFiles/Solvers.C headerFiles/Solvers.h headerFiles/Variables.h 
	${CXX} -c $< -o $@ ${CXXFLAGS}

objectFiles/RigidBodies.o : cppFiles/RigidBodies.C headerFiles/RigidBodies.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

objectFiles/LevelSet.o : cppFiles/LevelSet.C headerFiles/LevelSet.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

objectFiles/EOS.o : cppFiles/EOS.C headerFiles/EOS.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

objectFiles/eulerTests.o : cppFiles/eulerTests.C headerFiles/eulerTests.h 
	${CXX} -c $< -o $@ ${CXXFLAGS}

objectFiles/Variables.o : cppFiles/Variables.C headerFiles/Variables.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

objectFiles/Particles.o : cppFiles/Particles.C headerfiles/Particles.h
	${CXX} -c $< -o $@ ${CXXFLAGS}

main.o : main.C
	${CXX} -c $< -o $@ ${CXXFLAGS}

main : main.o   objectFiles/Solvers.o objectFiles/EOS.o objectFiles/eulerTests.o objectFiles/LevelSet.o objectFiles/Variables.o objectFiles/RigidBodies.o objectFiles/Particles.o
	g++ $^ -lm -o $@  

memcheck : main.o   objectFiles/Solvers.o objectFiles/EOS.o objectFiles/eulerTests.o objectFiles/LevelSet.o objectFiles/Variables.o objectFiles/RigidBodies.o objectFiles/Particles.o
	g++ -g  $^ -o $@



