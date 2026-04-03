# Makefile for Program TJ

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""
 
BUILD_DIR       = ../bin
LIB_DIR         = ../lib
UTILITY_DIR     = ../Utility
EQUILIBRIUM_DIR = ../Equilibrium
LAYER_DIR       = ../Layer
ISLAND_DIR      = ../Island

EQUILIBRIUM_LIB = $(LIB_DIR)/libEquilibrium.a
LAYER_LIB       = $(LIB_DIR)/libLayer.a
ISLAND_LIB      = $(LIB_DIR)/libIsland.a

C_HEADER = TJ.h
M_SOURCE = main.cpp
C_SOURCE = TJ.cpp ReadEquilibrium.cpp Interpolate.cpp Rational.cpp Matrix.cpp Vacuum.cpp Armadillo.cpp\
	ODESolve.cpp Resonant.cpp Dispersion.cpp Toroidal.cpp Ideal.cpp Netcdf.cpp RMP.cpp Visualize.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)
M_OBJECT = $(M_SOURCE:%.cpp=%.o)
U_OBJECT = ../Utility/Utility.o

LIBRARY    = ../lib/libTj.a
EXECUTABLE = ../bin/tj

all: 			$(LIBRARY) $(EXECUTABLE)

$(LIBRARY):	 	$(C_OBJECT) $(BUILD_DIR) $(EQUILIBRIUM_LIB) $(LAYER_LIB) $(ISLAND_LIB)
			$(C_ARCHIVER) $(LIBRARY) $(C_OBJECT) $(U_OBJECT) 

$(EXECUTABLE):	 	$(M_OBJECT) $(C_OBJECT) $(BUILD_DIR) $(EQUILIBRIUM_LIB) $(LAYER_LIB) $(ISLAND_LIB)
			$(C_COMPILER) -o $(EXECUTABLE) $(M_OBJECT) $(C_OBJECT) $(U_OBJECT) $(C_LIBS) -L../lib -lEquilibrium -lLayer -lIsland -lEce

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $*.cpp

$(M_OBJECT):		$(M_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $(M_SOURCE)

$(U_OBJECT):
			$(MAKE) -C $(UTILITY_DIR)

$(EQUILIBRIUM_LIB):
			$(MAKE) -C $(EQUILIBRIUM_DIR)

$(LAYER_LIB):
			$(MAKE) -C $(LAYER_DIR)

$(ISLAND_LIB):
			$(MAKE) -C $(ISLAND_DIR)

clean: 
			rm -f $(M_OBJECT) $(C_OBJECT) $(EXECUTABLE) $(LIBRARY)

