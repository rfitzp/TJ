# Makefile for library Layer

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

BUILD_DIR = ../lib

C_HEADER = Layer.h 
C_SOURCE = Layer.cpp 
C_OBJECT = $(C_SOURCE:%.cpp=%.o)

EXECUTABLE = ../lib/libLayer.a

$(EXECUTABLE):	 	$(C_OBJECT)  $(BUILD_DIR)
			$(C_ARCHIVER) $(EXECUTABLE) $(C_OBJECT)

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) $*.cpp

clean: 
			rm -f $(C_OBJECT) $(EXECUTABLE)

