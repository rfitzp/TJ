# Makefile for library Equilibrium

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

BUILD_DIR   = ../lib
LIB_DIR     = ../lib
UTILITY_DIR = ../Utility

C_HEADER = LightEquilibrium.h Equilibrium.h 
C_SOURCE = LightEquilibrium.cpp Equilibrium.cpp Netcdf.cpp EFIT.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o) 
U_OBJECT = ../Utility/Utility.o ../Utility/Faddeeva.o

EXECUTABLE = ../lib/libEquilibrium.a

$(EXECUTABLE):	 	$(C_OBJECT) $(BUILD_DIR)
			$(C_ARCHIVER) $(EXECUTABLE) $(C_OBJECT) $(U_OBJECT)

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $*.cpp

$(U_OBJECT):
			$(MAKE) -C $(UTILITY_DIR)

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

