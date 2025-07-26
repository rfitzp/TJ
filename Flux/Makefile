# Makefile for Program FLUX

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

BUILD_DIR   = ../bin
UTILITY_DIR = ../Utility

F_SOURCE = gFileRead.f90
F_OBJECT = $(F_SOURCE:%.f90=%.o)

C_HEADER = Flux.h
C_SOURCE = main.cpp Flux.cpp Stage2.cpp Interpolate.cpp Integrate.cpp Plasma.cpp 
C_OBJECT = $(C_SOURCE:%.cpp=%.o)
U_OBJECT = ../Utility/Utility.o ../Utility/Faddeeva.o

EXECUTABLE = ../bin/flux

$(EXECUTABLE):	 	$(C_OBJECT) $(F_OBJECT) $(BUILD_DIR)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(F_OBJECT) $(U_OBJECT) $(C_LIBS) $(F_LIBS)

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(F_OBJECT):		$(F_SOURCE)
			$(F_COMPILER) $(F_FLAGS) -o $*.o -c $(F_INCLUDE) $*.f90

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE)  -I../include $*.cpp

$(U_OBJECT):
			$(MAKE) -C $(UTILITY_DIR)

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)

