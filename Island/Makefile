# Makefile for program Island

include ../Machines.in

GIT_HASH     = `git rev-parse HEAD`
COMPILE_TIME = `date -u +'%Y-%m-%d %H:%M:%S UTC'`
GIT_BRANCH   = `git branch | grep "^\*" | sed 's/^..//'`
export VERSION_FLAGS=-DGIT_HASH="\"$(GIT_HASH)\"" -DCOMPILE_TIME="\"$(COMPILE_TIME)\"" -DGIT_BRANCH="\"$(GIT_BRANCH)\""

BUILD_DIR = ../lib
UTILITY_DIR = ../Utility

C_HEADER = Island.h 
C_SOURCE = Island.cpp
M_SOURCE = main.cpp
C_OBJECT = $(C_SOURCE:%.cpp=%.o)
M_OBJECT = $(M_SOURCE:%.cpp=%.o)
U_OBJECT = ../Utility/Utility.o

LIBRARY    = ../lib/libIsland.a
EXECUTABLE = ../bin/island

all: 			$(LIBRARY) $(EXECUTABLE)

$(LIBRARY):	 	$(C_OBJECT)  $(U_OBJECT) $(BUILD_DIR)
			$(C_ARCHIVER) $(LIBRARY) $(C_OBJECT) $(U_OBJECT)

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
			rm -f $(C_OBJECT) $(LIBRARY) $(EXECUTABLE)

