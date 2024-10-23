# Makefile for Program TEAR

include ../Machines.in

C_HEADER = Tear.h
C_SOURCE = main.cpp Tear.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../bin/tear

$(EXECUTABLE):	 	$(C_OBJECT) 
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS) 

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)

