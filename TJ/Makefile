# Makefile for Program TJ

include ../Machines.in

C_HEADER = TJ.h 
C_SOURCE = main.cpp TJ.cpp ReadEquilibrium.cpp Interpolate.cpp ZeroFind.cpp Rational.cpp\
	Matrix.cpp Vacuum.cpp Armadillo.cpp Integrate.cpp ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp Netcdf.cpp 
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../bin/tj

$(EXECUTABLE):	 	$(C_OBJECT) 
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS) -L../lib -lEquilibrium -lLayer

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

