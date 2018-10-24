cxx = g++
cxxflags = -Wall -Wextra


FDschemes.O : FDschemes.C FDschemes.H
	${cxx} -c FDschemes.C -o  FDschemes.O ${cxxflags}
FDmain.O : FDmain.C
	${cxx} -c FDmain.C -o FDmain.O ${cxxflags}
FDsolver : FDmain.O FDschemes.O
	${cxx} FDmain.O FDschemes.O -o FDsolver



