# Makefile for Program TJ

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""
 
BUILD_DIR       = ../bin
LIB_DIR         = ../lib
EQUILIBRIUM_DIR = ../Equilibrium
LAYER_DIR       = ../Layer

EQUILIBRIUM_LIB = $(LIB_DIR)/libEquilibrium.a
LAYER_LIB       = $(LIB_DIR)/libLayer.a

C_HEADER = TJ.h 
C_SOURCE = main.cpp TJ.cpp ReadEquilibrium.cpp Interpolate.cpp ZeroFind.cpp Rational.cpp\
	Matrix.cpp Vacuum.cpp Armadillo.cpp Integrate.cpp ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp\
	Ideal.cpp Netcdf.cpp 
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../bin/tj

$(EXECUTABLE):	 	$(C_OBJECT) $(BUILD_DIR) $(EQUILIBRIUM_LIB) $(LAYER_LIB)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS) -L../lib -lEquilibrium -lLayer

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(EQUILIBRIUM_LIB):
			$(MAKE) -C $(EQUILIBRIUM_DIR)

$(LAYER_LIB):
			$(MAKE) -C $(LAYER_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

