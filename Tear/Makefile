# Makefile for Program TEAR

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

BUILD_DIR = ../bin

C_HEADER = Tear.h
C_SOURCE = main.cpp Tear.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../bin/tear

$(EXECUTABLE):	 	$(C_OBJECT) $(BUILD_DIR)
			$(C_COMPILER) -o $(EXECUTABLE) $(C_OBJECT) $(C_LIBS)

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(F_OBJECT) $(C_OBJECT) $(EXECUTABLE)

