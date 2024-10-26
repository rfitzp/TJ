# Makefile for library Layer

include ../Machines.in

C_HEADER = Layer.h 
C_SOURCE = Layer.cpp 
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../lib/libLayer.a

$(EXECUTABLE):	 	$(C_OBJECT) 
			$(C_ARCHIVER) $(EXECUTABLE) $(C_OBJECT) 

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

