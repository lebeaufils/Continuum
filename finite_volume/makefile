cxx = g++
cxxflags = -Wall -Wextra


FVschemes.O : FVschemes.C FVschemes.H
	${cxx} -c FVschemes.C -o  FVschemes.O ${cxxflags}
FVmain.O : FVmain.C
	${cxx} -c FVmain.C -o FVmain.O ${cxxflags}
FVsolver : FVmain.O FVschemes.O
	${cxx} FVmain.O FVschemes.O -o FVsolver



