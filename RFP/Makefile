# Makefile for Program NONLINEAR

C_COMPILER = clang++ -O3
C_FLAGS    = -w 
C_INCLUDE  = -I/opt/local/include
C_LIBS     = -L/opt/local/lib -L/opt/local/lib -L/opt/local/lib/libgcc -lgfortran.5 -lgsl -lgslcblas

C_HEADER = NonLinear.h
C_OBJECT = main.o NonLinear.o
C_SOURCE = main.cpp NonLinear.cpp

EXECUTABLE = nonlinear

$(EXECUTABLE):	 	$(C_OBJECT)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

