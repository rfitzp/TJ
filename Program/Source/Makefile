# Makefile for Program TJ

C_COMPILER = clang++ 
C_FLAGS    = -w -O3 -std=c++11 
C_INCLUDE  = -I/opt/local/include
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgsl -lgslcblas -lnetcdf_c++4 -larmadillo -llapack

C_HEADER = TJ.h 
C_SOURCE = main.cpp TJ.cpp Equilibrium.cpp ReadEquilibrium.cpp Interpolate.cpp ZeroFind.cpp Rational.cpp Matrix.cpp\
	Vacuum.cpp Armadillo.cpp Integrate.cpp ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp Netcdf.cpp Layer.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../tj

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

