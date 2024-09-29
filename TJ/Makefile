# Makefile for Program TJ

include ../Machines.in

C_HEADER = LightEquilibrium.h Equilibrium.h TJ.h Layer.h
C_SOURCE = main.cpp TJ.cpp LightEquilibrium.cpp Equilibrium.cpp ReadEquilibrium.cpp Interpolate.cpp ZeroFind.cpp Rational.cpp\
	Matrix.cpp Vacuum.cpp Armadillo.cpp Integrate.cpp ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp Netcdf.cpp Layer.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../tj_

$(EXECUTABLE):	 	$(C_OBJECT) 
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

