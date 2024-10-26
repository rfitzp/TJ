# Makefile for library Equilibrium

include ../Machines.in

C_HEADER = LightEquilibrium.h Equilibrium.h 
C_SOURCE = LightEquilibrium.cpp Equilibrium.cpp Netcdf.cpp EFIT.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../lib/libEquilibrium.a

$(EXECUTABLE):	 	$(C_OBJECT) 
			$(C_ARCHIVER) $(EXECUTABLE) $(C_OBJECT) 

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

