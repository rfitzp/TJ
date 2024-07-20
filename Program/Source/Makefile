# Makefile for Program TJ

F_COMPILER = gfortran
F_FLAGS    =
F_INCLUDE  = -I/opt/local/include

C_COMPILER = clang++ 
C_FLAGS    = -w -O3 -std=c++11 
C_INCLUDE  = -I/opt/local/include
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgfortran.5 -lgsl -lgslcblas -lnetcdf_c++4 -larmadillo -llapack

F_SOURCE = NameListEquilibrium.f90 NameListTJ.f90
F_OBJECT = $(F_SOURCE:%.f90=%.o)

C_HEADER = TJ.h 
C_SOURCE = main.cpp TJ.cpp Equilibrium.cpp ReadEquilibrium.cpp Interpolate.cpp ZeroFind.cpp Rational.cpp Matrix.cpp\
	Vacuum.cpp Armadillo.cpp Integrate.cpp ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp Netcdf.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../tj

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(F_OBJECT) $(EXECUTABLE)

