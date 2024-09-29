# Makefile for Program FLUX

include ../Machines.in

F_SOURCE = gFileRead.f90
F_OBJECT = $(F_SOURCE:%.f90=%.o)

C_HEADER = Flux.h
C_SOURCE = main.cpp Flux.cpp Stage2.cpp Interpolate.cpp Integrate.cpp Plasma.cpp 
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../flux_

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(C_LIBS) $(F_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)

