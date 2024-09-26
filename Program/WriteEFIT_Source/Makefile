# Makefile for Program WriteEFIT

include ../Machines.in

F_SOURCE = WriteEFIT.f90
F_OBJECT = $(F_SOURCE:%.f90=%.o)

EXECUTABLE = ../wefit

$(EXECUTABLE):	 	$(F_OBJECT) 
			$(F_COMPILER) -o $(EXECUTABLE) $(F_OBJECT) $(F_LIBS)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

clean: 
			rm -f $(F_OBJECT) $(EXECUTABLE)

