# Makefile for Program TJ

F_COMPILER = gfortran
F_FLAGS    =
F_INCLUDE  = -I/opt/local/include

C_COMPILER = clang++ -O3 -std=c++11 
C_FLAGS    = -w 
C_INCLUDE  = -I/opt/local/include
C_LIBS     = -L/opt/local/lib -L/opt/local/lib/libgcc -lgfortran.5 -lgsl -lgslcblas -lnetcdf_c++4 -larmadillo -llapack

F_OBJECT = NameListEquilibrium.o NameListTJ.o 
F_SOURCE = NameListEquilibrium.f90 NameListTJ.f90

C_HEADER = TJ.h 
C_OBJECT = main.o TJ.o Equilibrium.o ReadEquilibrium.o Interpolate.o ZeroFind.o Rational.o Matrix.o\
	Vacuum.o Armadillo.o Integrate.o ODESolve.o Resonant.o Dispersion.o Toroidal.o
C_SOURCE = main.cpp TJ.cpp Equilibrium.cpp ReadEquilibrium.cpp Interpolate.cpp ZeroFind.cpp Rational.cpp Matrix.cpp\
	Vacuum.cpp Armadillo.cpp Integrate.cpp ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp

EXECUTABLE = ../tj

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(C_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(F_OBJECT) $(EXECUTABLE)

