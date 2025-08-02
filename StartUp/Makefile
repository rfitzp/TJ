# Makefile for Program StartUp

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""
 
BUILD_DIR       = ../bin
UTILITY_DIR     = ../Utility

C_HEADER = StartUp.h 
C_SOURCE = StartUp.cpp
M_SOURCE = main.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)
M_OBJECT = $(M_SOURCE:%.cpp=%.o)
U_OBJECT = ../Utility/Utility.o ../Utility/Faddeeva.o

EXECUTABLE = ../bin/startup

all: 			$(EXECUTABLE)

$(EXECUTABLE):	 	$(M_OBJECT) $(C_OBJECT) $(U_OBJECT) $(BUILD_DIR) 
			$(C_COMPILER) -o $(EXECUTABLE) $(M_OBJECT) $(C_OBJECT) $(U_OBJECT) $(C_LIBS) 

$(BUILD_DIR):
			@mkdir -p $(BUILD_DIR)

$(C_OBJECT):		$(C_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $*.cpp

$(M_OBJECT):		$(M_SOURCE) $(C_HEADER)
			$(C_COMPILER) $(C_FLAGS) -o $*.o -c $(C_INCLUDE) -I../include $(M_SOURCE)

$(U_OBJECT):
			$(MAKE) -C $(UTILITY_DIR)

clean: 
			rm -f $(C_OBJECT)  $(EXECUTABLE)

