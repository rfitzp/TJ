# Makefile for Program Vertical

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""
 
BUILD_DIR       = ../bin
LIB_DIR         = ../lib
UTILITY_DIR     = ../Utility
EQUILIBRIUM_DIR = ../Equilibrium

EQUILIBRIUM_LIB = $(LIB_DIR)/libEquilibrium.a

C_HEADER = Vertical.h 
C_SOURCE = main.cpp Vertical.cpp ReadEquilibrium.cpp Interpolate.cpp Netcdf.cpp Toroidal.cpp Vacuum.cpp\
		Armadillo.cpp Matrix.cpp ODESolve.cpp Ideal.cpp Visualize.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)
U_OBJECT = ../Utility/Utility.o

EXECUTABLE = ../bin/vertical

$(EXECUTABLE):	 	$(C_OBJECT) $(BUILD_DIR) $(EQUILIBRIUM_LIB) 
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(U_OBJECT) $(C_LIBS) -L../lib -lEquilibrium

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $*.cpp

$(U_OBJECT):
			$(MAKE) -C $(UTILITY_DIR)

$(EQUILIBRIUM_LIB):
			$(MAKE) -C $(EQUILIBRIUM_DIR)

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

